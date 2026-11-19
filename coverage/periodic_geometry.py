"""Exact rational coverage geometry on a cylinder of period P.

The infinite periodic plane set is handled with a finite, symmetric tiling.
Every input region is shifted by an integer number of periods so that its
minimum x lies in the fundamental strip [0, P); five tiles per region then
cover [-2P, 4P), leaving a full period of margin around the working zone
[-P, 2P].  Boundary walks are traced on the quotient cylinder (vertices are
identified modulo (P, 0)) and lifted back to the plane while tracing, so the
integer winding number falls out of the accumulated displacement.

No artificial cut edge is ever created: the seam only appears as a modular
vertex identity, never as a segment of its own.
"""
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key

from .arrangement import direction_compare, face_witness
from .rings import region_contains
from .segments import AtomicSegment, SourceSpan, split_segments

TILE_OFFSETS = (-2, -1, 0, 1, 2)
EDGE_LIMIT = 80


@dataclass(frozen=True)
class Tile:
    region_index: int
    shift: int
    bbox: tuple
    rings: tuple


@dataclass(frozen=True)
class WalkSegment:
    delta: tuple
    sources: tuple


@dataclass(frozen=True)
class PeriodicWalk:
    boundary_id: str
    start: tuple
    segments: tuple
    winding: int


@dataclass(frozen=True)
class PeriodicComponent:
    component_id: str
    walks: tuple


@dataclass(frozen=True)
class PeriodicGeometry:
    period: Fraction
    components: tuple


class _UnionFind:
    def __init__(self, size):
        self._parent = list(range(size))

    def root(self, i):
        while self._parent[i] != i:
            self._parent[i] = self._parent[self._parent[i]]
            i = self._parent[i]
        return i

    def join(self, values):
        values = tuple(values)
        if len(values) < 2:
            return
        first = self.root(values[0])
        for value in values[1:]:
            self._parent[self.root(value)] = first


def canonical_shift(rings, period):
    """Integer k so that min(x) + kP lies in [0, P)."""
    min_x = min(point[0] for ring in rings for point in ring)
    return -(min_x // period)


def validate_periodic_inputs(revisions, period):
    for revision in revisions:
        xs = [point[0] for ring in revision.rings for point in ring]
        if max(xs) - min(xs) >= period:
            raise ValueError('region horizontal span must be strictly less than period')
    edge_count = sum(len(ring) for revision in revisions for ring in revision.rings)
    if edge_count > EDGE_LIMIT:
        raise ValueError('periodic batch exceeds the fixed edge limit')


def _shift_rings(rings, shift, period):
    offset = shift * period
    return tuple(tuple((point[0] + offset, point[1]) for point in ring) for ring in rings)


def _bbox(rings):
    xs = [point[0] for ring in rings for point in ring]
    ys = [point[1] for ring in rings for point in ring]
    return min(xs), max(xs), min(ys), max(ys)


def build_tiles(revisions, period):
    tiles = []
    shifts = []
    for index, revision in enumerate(revisions):
        base = canonical_shift(revision.rings, period)
        shifts.append(base)
        for offset in TILE_OFFSETS:
            shift = base + offset
            rings = _shift_rings(revision.rings, shift, period)
            tiles.append(Tile(index, shift, _bbox(rings), rings))
    return tuple(tiles), tuple(shifts)


def tiled_source_segments(revisions, tiles):
    segments = []
    for tile in tiles:
        revision = revisions[tile.region_index]
        for ring_index, ring in enumerate(tile.rings):
            for edge_index, (a, b) in enumerate(zip(ring, ring[1:] + ring[:1])):
                source = SourceSpan(revision.region_id, revision.revision,
                                    ring_index, edge_index, Fraction(0), Fraction(1))
                segments.append(AtomicSegment(a, b, (source,)))
    return tuple(segments)


def covering_tiles(witness, tiles):
    wx, wy = witness
    result = []
    for index, tile in enumerate(tiles):
        min_x, max_x, min_y, max_y = tile.bbox
        if min_x <= wx <= max_x and min_y <= wy <= max_y and region_contains(witness, tile.rings):
            result.append(index)
    return tuple(result)


def classify(segments, tiles):
    sides = []
    for segment in segments:
        pair = []
        for side in (1, -1):
            witness = face_witness(segment, side, segments)
            pair.append(covering_tiles(witness, tiles))
        sides.append(tuple(pair))
    return tuple(sides)


def mod_point(point, period):
    return point[0] - period * (point[0] // period), point[1]


def _source_key(source):
    return (source.region_id, source.revision, source.ring_index,
            source.edge_index, source.t0, source.t1)


def build_periodic_geometry(revisions, period):
    validate_periodic_inputs(revisions, period)
    revisions = tuple(sorted(revisions, key=lambda r: (r.region_id, r.revision)))
    tiles, _ = build_tiles(revisions, period)
    atomic = split_segments(tiled_source_segments(revisions, tiles))
    sides = classify(atomic, tiles)

    finder = _UnionFind(len(revisions))
    for left, right in sides:
        left_regions = {tiles[i].region_index for i in left}
        right_regions = {tiles[i].region_index for i in right}
        finder.join(left_regions)
        finder.join(right_regions)
        if left_regions and right_regions:
            finder.join(left_regions | right_regions)

    # Exposed atomic segments.  A segment is a true cylinder boundary only if
    # its coverage classification is reliable, which the five-tile window
    # guarantees inside [-P, 2P]; every planar lift of such an edge is kept
    # for continuous lift tracing across the seam.
    boundary_keys = set()
    for segment, (left, right) in zip(atomic, sides):
        covered_left = bool(left)
        if covered_left == bool(right):
            continue
        oriented = segment if covered_left else segment.reversed()
        if all(-period <= point[0] <= 2 * period for point in (oriented.start, oriented.end)):
            boundary_keys.add((mod_point(oriented.start, period),
                               oriented.end[0] - oriented.start[0],
                               oriented.end[1] - oriented.start[1]))

    grouped = {}
    for segment, (left, right) in zip(atomic, sides):
        covered_left = bool(left)
        if covered_left == bool(right):
            continue
        oriented = segment if covered_left else segment.reversed()
        key = (mod_point(oriented.start, period), oriented.end[0] - oriented.start[0],
               oriented.end[1] - oriented.start[1])
        if key not in boundary_keys:
            continue
        covering = left if covered_left else right
        root = finder.root(tiles[covering[0]].region_index)
        entry = grouped.setdefault(key, {'sources': set(), 'roots': set(), 'lifts': {}})
        entry['sources'].update(oriented.sources)
        entry['roots'].add(root)
        level_fraction = (oriented.start[0] - key[0][0]) / period
        if level_fraction.denominator != 1 or oriented.start[1] != key[0][1]:
            raise ValueError('inconsistent periodic edge lift')
        level = int(level_fraction)
        if level in entry['lifts']:
            raise ValueError('inconsistent periodic edge lift')
        entry['lifts'][level] = oriented

    edges = []
    for index, (key, entry) in enumerate(sorted(grouped.items())):
        start, dx, dy = key
        roots = tuple(sorted(entry['roots']))
        if len({finder.root(r) for r in roots}) != 1:
            raise ValueError('periodic edge joins distinct coverage components')
        edges.append({
            'index': index,
            'start': start,
            'delta': (dx, dy),
            'sources': tuple(sorted(entry['sources'], key=_source_key)),
            'root': roots[0],
            'lifts': dict(entry['lifts']),
        })

    cycles = trace_cycles(edges, period)
    walks = []
    for cycle in cycles:
        walk = lift_cycle(cycle, period)
        walks.append((cycle[0]['root'], walk))

    component_roots = {}
    for root, walk in walks:
        component_roots.setdefault(root, []).append(walk)
    if set(component_roots) != {finder.root(i) for i in range(len(revisions))}:
        raise ValueError('a coverage component has no boundary')

    component_keys = []
    for root, loop_walks in component_roots.items():
        loop_keys = tuple(sorted(walk_sort_key(walk) for walk in loop_walks))
        component_keys.append((loop_keys, root, loop_walks))
    component_keys.sort(key=lambda item: item[0])

    components = []
    boundary_number = 0
    for component_index, (_, root, loop_walks) in enumerate(component_keys):
        ordered = sorted(loop_walks, key=walk_sort_key)
        walk_objects = []
        for walk in ordered:
            walk_objects.append(PeriodicWalk(
                'b' + str(boundary_number), walk.start, walk.segments, walk.winding))
            boundary_number += 1
        components.append(PeriodicComponent('c' + str(component_index), tuple(walk_objects)))
    return PeriodicGeometry(period, tuple(components))


def trace_cycles(edges, period):
    fan = {}
    for edge in edges:
        fan.setdefault(edge['start'], []).append(edge)
    for outgoing in fan.values():
        directions = [edge['delta'] for edge in outgoing]
        if len(set(directions)) != len(directions):
            raise ValueError('ambiguous periodic vertex fan')

    successors = {}
    for edge in edges:
        end = mod_point((edge['start'][0] + edge['delta'][0],
                         edge['start'][1] + edge['delta'][1]), period)
        vectors = [(other['delta'], other['index']) for other in fan[end]]
        reverse = (-edge['delta'][0], -edge['delta'][1])
        vectors.append((reverse, -1))
        vectors.sort(key=cmp_to_key(lambda a, b: direction_compare(a[0], b[0])))
        position = next(k for k, item in enumerate(vectors) if item[1] == -1)
        successors[edge['index']] = vectors[position - 1][1]
    if sorted(successors.values()) != list(range(len(edges))):
        raise ValueError('periodic boundary continuation is not a permutation')

    by_index = {edge['index']: edge for edge in edges}
    remaining = set(by_index)
    cycles = []
    while remaining:
        start_index = min(remaining)
        current = start_index
        cycle = []
        while current in remaining:
            remaining.remove(current)
            cycle.append(by_index[current])
            current = successors[current]
        if current != start_index:
            raise ValueError('periodic boundary joins another walk')
        cycles.append(tuple(cycle))
    return tuple(cycles)


def lift_cycle(cycle, period):
    first = cycle[0]
    start = first['start']
    if 0 not in first['lifts']:
        raise ValueError('periodic boundary leaves the tiling window')
    point = start
    level = 0
    lifted = []
    for position, edge in enumerate(cycle):
        lift = edge['lifts'].get(level)
        if lift is None or lift.start != point:
            raise ValueError('periodic boundary leaves the tiling window')
        next_start = cycle[(position + 1) % len(cycle)]['start']
        if mod_point(lift.end, period) != next_start:
            raise ValueError('periodic boundary lift mismatch')
        lifted.append(WalkSegment(edge['delta'], edge['sources']))
        # The planar endpoint is the continuous lift.  Crossing the seam keeps
        # the planar point while the next edge is indexed by its lift level.
        point = lift.end
        if position + 1 < len(cycle):
            next_level_fraction = (point[0] - next_start[0]) / period
            if next_level_fraction.denominator != 1 or point[1] != next_start[1]:
                raise ValueError('periodic boundary lift mismatch')
            level = int(next_level_fraction)
    if point[1] != start[1]:
        raise ValueError('periodic walk changes y level')
    span = point[0] - start[0]
    winding_fraction = span / period
    if winding_fraction.denominator != 1:
        raise ValueError('non-integer periodic winding')
    winding = int(winding_fraction)
    rotation = canonical_rotation(cycle)
    rotated = cycle[rotation:] + cycle[:rotation]
    segments = tuple(WalkSegment(edge['delta'], edge['sources']) for edge in rotated)
    return PeriodicWalk('', rotated[0]['start'], segments, winding)


def canonical_rotation(cycle):
    def key(rotation):
        return tuple((edge['start'], edge['delta'], edge['sources']) for edge in
                     cycle[rotation:] + cycle[:rotation])
    return min(range(len(cycle)), key=key)


def walk_sort_key(walk):
    return tuple((walk.start, segment.delta, segment.sources) for segment in walk.segments)


def walk_area(walk):
    # Green's theorem on the cylinder with the globally defined form -y dx:
    # d(-y dx) = dx ^ dy.  Valid for every winding, including winding +/-1.
    total = Fraction(0)
    y = walk.start[1]
    for segment in walk.segments:
        dx, dy = segment.delta
        next_y = y + dy
        total -= (y + next_y) * dx / 2
        y = next_y
    return total


def component_area(component):
    return sum((walk_area(walk) for walk in component.walks), Fraction())


def periodic_area_items(geometry):
    items = [{'component_id': component.component_id, 'area': component_area(component)}
             for component in geometry.components]
    if any(item['area'] <= 0 for item in items):
        raise ValueError('nonpositive periodic component area')
    return items


def periodic_total_area(geometry):
    return sum((component_area(component) for component in geometry.components), Fraction())


def segment_length_squared(segment):
    return segment.delta[0] ** 2 + segment.delta[1] ** 2


def periodic_boundary_measurements(geometry):
    result = []
    for component in geometry.components:
        for walk in component.walks:
            result.append({'boundary_id': walk.boundary_id,
                           'length_squared': [segment_length_squared(segment) for segment in walk.segments]})
    return result
