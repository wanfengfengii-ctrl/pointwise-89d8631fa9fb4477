from fractions import Fraction
from dataclasses import asdict
from .reading import read_report
from .boundary_measure import find_boundary_segment
from .intersections import interpolate
from .periodic_codec import periodic_geometry_from_data

def validate_sample_parameter(t):
    if type(t) is not Fraction or not 0 <= t <= 1:
        raise ValueError('t must be a Fraction in [0,1]')
    return t

def interpolate_source_span(source, t):
    result = asdict(source)
    result['parameter'] = source.t0 + t*(source.t1-source.t0)
    return result

def find_periodic_segment(geometry, boundary_id, segment_index):
    if type(segment_index) is not int or isinstance(segment_index, bool) or segment_index < 0:
        raise KeyError((boundary_id, segment_index))
    for component in geometry.components:
        for walk in component.walks:
            if walk.boundary_id == boundary_id:
                if segment_index >= len(walk.segments):
                    raise KeyError((boundary_id, segment_index))
                return walk, walk.segments[segment_index]
    raise KeyError(boundary_id)

def sample_boundary(report_id, boundary_id, segment_index, t):
    validate_sample_parameter(t)
    report = read_report(report_id)
    if report['format'] == 2:
        geometry = periodic_geometry_from_data(report['components'], report['period'])
        walk, segment = find_periodic_segment(geometry, boundary_id, segment_index)
        lifted_start = walk.start
        for earlier in walk.segments[:segment_index]:
            lifted_start = (lifted_start[0] + earlier.delta[0],
                            lifted_start[1] + earlier.delta[1])
        lifted_end = (lifted_start[0] + segment.delta[0],
                      lifted_start[1] + segment.delta[1])
        point = interpolate(lifted_start, lifted_end, t)
        period = report['period']
        point = [point[0] - period * (point[0] // period), point[1]]
        return {'point': point,
                'sources': [interpolate_source_span(source, t) for source in segment.sources]}
    from .coverage_union import geometry_from_data
    geometry = geometry_from_data(report['components'])
    segment = find_boundary_segment(geometry, boundary_id, segment_index)
    return {'point': list(interpolate(segment.start, segment.end, t)),
            'sources': [interpolate_source_span(source, t) for source in segment.sources]}
