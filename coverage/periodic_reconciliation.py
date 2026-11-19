"""Fixed-consumption verification for format 2 periodic reports.

Only the saved representation is consumed: the cylinder coverage itself is
never recomputed.  Checks cover period closure and integer winding, exact
recomputed areas and per-segment squared lengths, and one-to-one periodic
correspondence between every saved source interval and the original edges in
the saved inputs, including all periodic coincident contributors.
"""
from fractions import Fraction

from .intersections import parameter_on_segment, point_on_segment, segment_intersections
from .rings import ring_edges
from .periodic_geometry import (periodic_area_items,
                                periodic_boundary_measurements,
                                periodic_total_area, segment_length_squared)


def verify_walk_closure(geometry):
    for component in geometry.components:
        for walk in component.walks:
            point = walk.start
            for segment in walk.segments:
                point = (point[0] + segment.delta[0], point[1] + segment.delta[1])
            if point[1] != walk.start[1]:
                raise ValueError('periodic walk does not return to its y level')
            span = point[0] - walk.start[0]
            winding_fraction = span / geometry.period
            if winding_fraction.denominator != 1 or int(winding_fraction) != walk.winding:
                raise ValueError('wrong periodic winding')


def verify_measurements(report, geometry):
    expected = {
        'areas': periodic_area_items(geometry),
        'total_area': periodic_total_area(geometry),
        'boundary': periodic_boundary_measurements(geometry),
    }
    if report['measurements'] != expected:
        raise ValueError('periodic measurement reconciliation failed')
    areas = {item['component_id']: item['area'] for item in expected['areas']}
    if sum(areas.values(), Fraction()) != expected['total_area']:
        raise ValueError('periodic total area mismatch')
    saved = {item['boundary_id']: item['length_squared']
             for item in report['measurements']['boundary']}
    for component in geometry.components:
        for walk in component.walks:
            expected_squares = [segment_length_squared(segment) for segment in walk.segments]
            if saved.get(walk.boundary_id) != expected_squares:
                raise ValueError('periodic squared length mismatch')


def verify_periodic_sources(report, geometry):
    period = geometry.period
    originals = []
    for revision in report['inputs']:
        for ring_index, ring in enumerate(revision['rings']):
            for edge_index, (a, b) in enumerate(ring_edges(tuple(tuple(p) for p in ring))):
                originals.append((revision['region_id'], revision['revision'],
                                  ring_index, edge_index, tuple(a), tuple(b)))
    # A lifted segment starts inside one period of [0, P); enumerate every
    # integer lift of each original edge that could reach that coordinate band.
    edge_min = min(a[0] for _, _, _, _, a, _ in originals)
    edge_max = max(b[0] for _, _, _, _, _, b in originals)
    low = int((-period - edge_min) // period) - 1
    high = int((2 * period - edge_max) // period) + 1
    shifts = tuple(range(low, high + 1))

    def lift_edge(a, b, shift):
        return (a[0] + shift * period, a[1]), (b[0] + shift * period, b[1])

    for component in geometry.components:
        for walk in component.walks:
            point = walk.start
            for segment in walk.segments:
                start, end = point, (point[0] + segment.delta[0], point[1] + segment.delta[1])
                expected = []
                for region_id, revision, ring_index, edge_index, a, b in originals:
                    matches = []
                    for shift in shifts:
                        la, lb = lift_edge(a, b, shift)
                        hits = segment_intersections(start, end, la, lb)
                        if len(hits) == 2:
                            if not (point_on_segment(start, la, lb) and point_on_segment(end, la, lb)):
                                raise ValueError('periodic boundary omits an original endpoint split')
                            matches.append((shift, parameter_on_segment(start, la, lb),
                                            parameter_on_segment(end, la, lb)))
                        elif hits and 0 < hits[0].first_parameter < 1:
                            raise ValueError('periodic boundary omits an intersection split')
                    if len(matches) > 1:
                        raise ValueError('an original edge matches multiple periodic lifts')
                    if matches:
                        shift, t0, t1 = matches[0]
                        la, lb = lift_edge(a, b, shift)
                        # The two interpolated source points must differ from
                        # the lifted segment endpoints by the same (kP, 0).
                        source_start = (la[0] + t0 * (lb[0] - la[0]),
                                        la[1] + t0 * (lb[1] - la[1]))
                        source_end = (la[0] + t1 * (lb[0] - la[0]),
                                      la[1] + t1 * (lb[1] - la[1]))
                        offset_x = start[0] - source_start[0]
                        if start[1] != source_start[1] or end[1] != source_end[1]:
                            raise ValueError('source interval does not follow the boundary y')
                        if (offset_x / period).denominator != 1 or end[0] - source_end[0] != offset_x:
                            raise ValueError('source lift is not an integer number of periods')
                        expected.append((region_id, revision, ring_index, edge_index, t0, t1))
                actual = tuple((source.region_id, source.revision, source.ring_index,
                                source.edge_index, source.t0, source.t1)
                               for source in segment.sources)
                if tuple(sorted(expected)) != actual:
                    raise ValueError('incorrect or incomplete periodic boundary sources')
                point = end


def reconcile_periodic_report(report):
    from .periodic_codec import validate_periodic_report_structure
    geometry = validate_periodic_report_structure(report)
    verify_walk_closure(geometry)
    verify_measurements(report, geometry)
    verify_periodic_sources(report, geometry)
    boundaries = [walk for component in geometry.components for walk in component.walks]
    return {
        'component_count': len(geometry.components),
        'boundary_count': len(boundaries),
        'segment_count': sum(len(walk.segments) for walk in boundaries),
        'total_area': periodic_total_area(geometry),
        'period': geometry.period,
    }
