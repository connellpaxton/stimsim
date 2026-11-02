"""Fixed-step optimizations validated for the static published current-clamp protocol."""


def prune_zero_conductance_channels(sections):
    """Remove unused Nap_Et2/K_Tst only when all section conductances are zero.

    Do not use this when later code changes these conductances, records their
    gating states, or holds pointers into the removed mechanisms.
    Return the number of removed computational-segment instances per mechanism.
    """
    from neuron import h

    removed = {}
    for section in sections:
        for mechanism, parameter in [('Nap_Et2', 'gNap_Et2bar'), ('K_Tst', 'gK_Tstbar')]:
            if h.ismembrane(mechanism, sec=section) and all(
                getattr(segment, f'{parameter}_{mechanism}') == 0.
                for segment in section
            ):
                removed[mechanism] = removed.get(mechanism, 0) + section.nseg
                section.uninsert(mechanism)
    return removed
