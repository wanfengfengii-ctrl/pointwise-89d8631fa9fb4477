"""Structure validation and persistence geometry for format 2 reports."""
from fractions import Fraction

from .coverage_union import require_keys
from .periodic_geometry import (PeriodicComponent, PeriodicGeometry,
                                PeriodicWalk, WalkSegment)
from .region_validation import validate_region
from .rings import identifier
from .segments import SourceSpan

REPORT_FORMAT_VERSION = 2


def rational_pair(value):
    if not isinstance(value, list) or len(value) != 2 or any(type(v) is not Fraction for v in value):
        raise ValueError('invalid exact point or delta')
    return tuple(value)


def periodic_geometry_from_data(data, period):
    if not isinstance(data, list) or not data:
        raise ValueError('missing components')
    components, component_ids, boundary_ids = [], set(), set()
    for component in data:
        require_keys(component, ('component_id', 'boundaries'))
        identifier(component['component_id'])
        if component['component_id'] in component_ids or not isinstance(component['boundaries'], list) \
                or not component['boundaries']:
            raise ValueError('duplicate component or missing boundaries')
        component_ids.add(component['component_id'])
        walks = []
        for boundary in component['boundaries']:
            require_keys(boundary, ('boundary_id', 'start', 'segments', 'winding'))
            identifier(boundary['boundary_id'])
            if boundary['boundary_id'] in boundary_ids:
                raise ValueError('duplicate boundary identity')
            boundary_ids.add(boundary['boundary_id'])
            if type(boundary['winding']) is not int or isinstance(boundary['winding'], bool):
                raise ValueError('invalid winding')
            start = rational_pair(boundary['start'])
            if not 0 <= start[0] < period:
                raise ValueError('boundary start outside the fundamental strip')
            if not isinstance(boundary['segments'], list) or not boundary['segments']:
                raise ValueError('missing boundary segments')
            segments = []
            for segment in boundary['segments']:
                require_keys(segment, ('delta', 'sources'))
                delta = rational_pair(segment['delta'])
                if delta == (Fraction(0), Fraction(0)):
                    raise ValueError('zero boundary displacement')
                if not isinstance(segment['sources'], list) or not segment['sources']:
                    raise ValueError('missing boundary sources')
                sources = []
                for source in segment['sources']:
                    require_keys(source, ('region_id', 'revision', 'ring_index',
                                          'edge_index', 't0', 't1'))
                    identifier(source['region_id'])
                    for key, lower in (('revision', 1), ('ring_index', 0), ('edge_index', 0)):
                        if type(source[key]) is not int or source[key] < lower:
                            raise ValueError('invalid source identity')
                    if any(type(source[k]) is not Fraction or not 0 <= source[k] <= 1
                           for k in ('t0', 't1')) or source['t0'] == source['t1']:
                        raise ValueError('invalid source interval')
                    sources.append(SourceSpan(**source))
                if sources != sorted(set(sources)):
                    raise ValueError('duplicate or unordered sources')
                segments.append(WalkSegment(delta, tuple(sources)))
            walks.append(PeriodicWalk(boundary['boundary_id'], start, tuple(segments),
                                      boundary['winding']))
        components.append(PeriodicComponent(component['component_id'], tuple(walks)))
    return PeriodicGeometry(period, tuple(components))


def validate_periodic_report_structure(report):
    require_keys(report, ('format', 'report_id', 'batch_id', 'binding', 'inputs',
                          'period', 'components', 'measurements', 'receipt'))
    if type(report['format']) is not int or isinstance(report['format'], bool) \
            or report['format'] != REPORT_FORMAT_VERSION:
        raise ValueError('unsupported report format')
    identifier(report['report_id'])
    identifier(report['batch_id'])
    period = report['period']
    if type(period) is not Fraction or period <= 0:
        raise ValueError('period must be a positive Fraction')
    if not isinstance(report['binding'], list) or not report['binding'] \
            or not isinstance(report['inputs'], list):
        raise ValueError('missing fixed inputs')
    keys = []
    for member in report['binding']:
        require_keys(member, ('region_id', 'revision'))
        identifier(member['region_id'])
        if type(member['revision']) is not int or member['revision'] <= 0:
            raise ValueError('invalid revision')
        keys.append((member['region_id'], member['revision']))
    if keys != sorted(keys) or len({key[0] for key in keys}) != len(keys):
        raise ValueError('duplicate or unordered binding')
    snapshots = []
    for source in report['inputs']:
        require_keys(source, ('region_id', 'revision', 'rings'))
        identifier(source['region_id'])
        if type(source['revision']) is not int or source['revision'] <= 0:
            raise ValueError('invalid snapshot revision')
        rings = source['rings']
        if not isinstance(rings, list) or not rings:
            raise ValueError('missing input rings')
        for ring in rings:
            if not isinstance(ring, list):
                raise ValueError('invalid input ring')
            for point in ring:
                rational_pair(point)
        validated = validate_region(rings[0], rings[1:])
        if any(len(a) != len(b) for a, b in zip(validated.rings, rings)):
            raise ValueError('snapshot must omit closing duplicate')
        snapshots.append((source['region_id'], source['revision']))
    if snapshots != keys:
        raise ValueError('snapshot binding mismatch')
    geometry = periodic_geometry_from_data(report['components'], period)
    measurements = report['measurements']
    require_keys(measurements, ('areas', 'total_area', 'boundary'))
    if type(measurements['total_area']) is not Fraction \
            or not isinstance(measurements['areas'], list) \
            or not isinstance(measurements['boundary'], list):
        raise ValueError('missing measurement items')
    for item in measurements['areas']:
        require_keys(item, ('component_id', 'area'))
        if type(item['area']) is not Fraction or item['area'] <= 0:
            raise ValueError('invalid component area')
    for item in measurements['boundary']:
        require_keys(item, ('boundary_id', 'length_squared'))
        if not isinstance(item['length_squared'], list) or not item['length_squared'] \
                or any(type(v) is not Fraction for v in item['length_squared']):
            raise ValueError('invalid squared lengths')
    require_keys(report['receipt'], ('report_id', 'batch_id', 'status'))
    if report['receipt'] != {'report_id': report['report_id'],
                             'batch_id': report['batch_id'], 'status': 'published'}:
        raise ValueError('receipt binding mismatch')
    return geometry
