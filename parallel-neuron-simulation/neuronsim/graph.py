"""Count contacts and connected components of the observed directed graph."""

from collections import Counter


def graph_summary(nodes, edges):
    ids = {row['root_id'] for row in nodes}
    forward = {root: set() for root in ids}
    reverse = {root: set() for root in ids}
    outgoing = Counter()
    incoming = Counter()
    total = 0
    for row in edges:
        pre, post = (row['pre_pt_root_id'], row['post_pt_root_id'])
        count = int(row['n_syn'])
        if pre not in ids or post not in ids or count <= 0:
            raise ValueError('Invalid internal edge')
        forward[pre].add(post)
        reverse[post].add(pre)
        outgoing[pre] += count
        incoming[post] += count
        total += count
    visited = set()
    order = []
    for root in sorted(ids):
        if root in visited:
            continue
        visited.add(root)
        stack = [(root, iter(forward[root]))]
        while stack:
            node, neighbors = stack[-1]
            try:
                neighbor = next(neighbors)
            except StopIteration:
                order.append(node)
                stack.pop()
                continue
            if neighbor not in visited:
                visited.add(neighbor)
                stack.append((neighbor, iter(forward[neighbor])))

    def components(roots, neighbors):
        seen = set()
        sizes = []
        for root in roots:
            if root in seen:
                continue
            seen.add(root)
            todo = [root]
            size = 0
            while todo:
                node = todo.pop()
                size += 1
                for other in neighbors(node):
                    if other not in seen:
                        seen.add(other)
                        todo.append(other)
            sizes.append(size)
        return sorted(sizes, reverse=True)
    strong = components(reversed(order), lambda node: reverse[node])
    weak = components(sorted(ids), lambda node: forward[node] | reverse[node])
    return {'neurons': len(ids), 'anatomical_contacts': total, 'directed_pairs': len(edges), 'directed_pair_density': len(edges) / (len(ids) * (len(ids) - 1)) if len(ids) > 1 else 0, 'weak_component_sizes': weak, 'strong_component_sizes': strong, 'isolated_neurons': sorted((root for root in ids if not forward[root] and (not reverse[root]))), 'neurons_without_internal_input': sum((not reverse[root] for root in ids)), 'neurons_without_internal_output': sum((not forward[root] for root in ids)), 'incoming_contacts': dict(incoming), 'outgoing_contacts': dict(outgoing)}
