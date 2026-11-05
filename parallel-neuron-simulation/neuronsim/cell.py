"""Turn a recorded SWC tree into a NEURON cell with published biophysics."""

import itertools
from pathlib import Path
import re

import numpy as np

from .optimizations import prune_zero_conductance_channels
from .swc import load_swc


class PublishedSWCCell:
    """Keep the recorded branches; assign the supplied HOC regional profile."""

    _instances = itertools.count()

    def __init__(self, path, profile, *, original_discretization=False, coordinate_cache=True):
        from neuron import h

        self.path = Path(path)
        self.cell_id = self.path.stem
        self.instance = next(self._instances)
        rows, _, _ = load_swc(self.path)
        if not set(rows[:, 1]).issubset({1, 2, 3, 4}):
            raise ValueError(f"Unsupported SWC section types in {self.path}")
        self.swc_nodes = len(rows)

        # SIMULATION CORE: NEURON Import3D creates the connected cable sections.
        reader = h.Import3d_SWC_read()
        reader.input(str(self.path))
        importer = h.Import3d_GUI(reader, False)
        importer.instantiate(self)
        self.sections = list(self.all)
        if not getattr(self, "soma", None):
            raise ValueError(f"No soma in {self.path}")

        # These are the region names expected by the unchanged upstream HOC.
        for name, source in [
            ("all", None), ("somatic", "soma"), ("basal", "dend"),
            ("apical", "apic"), ("axonal", "axon"),
        ]:
            section_list = h.SectionList()
            selected = self.sections if source is None else list(getattr(self, source, []))
            for section in selected:
                section_list.append(sec=section)
            setattr(self, name, section_list)

        # SIMULATION CORE: nseg determines the electrical compartments solved.
        for section in self.sections:
            if original_discretization:
                # cACint209 uses Ra=100 and cm=1 in every region. These determine
                # lambda_f. Match the supplied script's nseg reduction factor 2.
                section.Ra = 100.0
                section.cm = 1.0
                full = int((section.L / (0.01 * h.lambda_f(100.0, sec=section)) + 0.999) / 2) * 2 + 1
                section.nseg = int(full / 9 / 2) * 2 + 1
            else:
                # Keep the earlier L23 adapter's discretization as its default.
                section.nseg = 1 + 2 * int(section.L / 40)

        # Set the origin before HOC evaluates distance-dependent channel density.
        self.soma[0].push()
        h.distance()
        h.pop_section()
        # SIMULATION CORE: insert channels and set region-specific conductances.
        profile.biophys(self)
        self.removed_zero_channels = prune_zero_conductance_channels(self.all)

        if not coordinate_cache:
            self.segments = [segment for section in self.sections for segment in section]
            self.segment_xyz_um = None
            return

        # Older contact-coordinate datasets use this nearest-segment lookup.
        self.segments = []
        coordinates = []
        for section in self.sections:
            points = np.array([
                [section.x3d(i), section.y3d(i), section.z3d(i)]
                for i in range(section.n3d())
            ])
            arcs = np.array([section.arc3d(i) for i in range(section.n3d())])
            if len(arcs) < 2 or arcs[-1] <= 0:
                raise ValueError(f"Invalid section geometry: {section.name()}")
            for segment in section:
                coordinates.append([
                    np.interp(segment.x * arcs[-1], arcs, points[:, axis])
                    for axis in range(3)
                ])
                self.segments.append(segment)
        self.segment_xyz_um = np.asarray(coordinates)

    def __str__(self):
        # Upstream distribute() executes HOC containing this object's identifier.
        identifier = re.sub(r"\W", "_", self.cell_id)
        return f"swc_{identifier}_{self.instance}"

    def nearest_segment(self, position_nm):
        position = np.asarray(position_nm, dtype=float) / 1000.0
        difference = self.segment_xyz_um - position
        distances_squared = np.einsum("ij,ij->i", difference, difference)
        index = int(distances_squared.argmin())
        return self.segments[index], float(np.sqrt(distances_squared[index]))
