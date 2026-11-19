from fractions import Fraction

from .batches import read_sealed_batch
from .database import current_archive
from .periodic_geometry import (build_periodic_geometry,
                                periodic_area_items, periodic_boundary_measurements,
                                periodic_total_area)
from .periodic_reconciliation import reconcile_periodic_report
from .publication import publication_kind, publish_report
from .rings import identifier


def normalize_period(period):
    if isinstance(period, bool) or not isinstance(period, (int, Fraction)):
        raise ValueError('period must be a positive integer or Fraction')
    value = Fraction(period)
    if value <= 0:
        raise ValueError('period must be strictly positive')
    return value


def geometry_data(geometry):
    return [{
        'component_id': component.component_id,
        'boundaries': [{
            'boundary_id': walk.boundary_id,
            'start': list(walk.start),
            'segments': [{
                'delta': list(segment.delta),
                'sources': [{
                    'region_id': source.region_id,
                    'revision': source.revision,
                    'ring_index': source.ring_index,
                    'edge_index': source.edge_index,
                    't0': source.t0,
                    't1': source.t1,
                } for source in segment.sources],
            } for segment in walk.segments],
            'winding': walk.winding,
        } for walk in component.walks],
    } for component in geometry.components]


def assemble_periodic_report(report_id, binding, period, geometry):
    report = {
        'format': 2,
        'report_id': report_id,
        'batch_id': binding.batch_id,
        'binding': [{'region_id': member.region_id, 'revision': member.revision}
                    for member in binding.members],
        'inputs': [{
            'region_id': revision.region_id,
            'revision': revision.revision,
            'rings': [[list(point) for point in ring] for ring in revision.rings],
        } for revision in binding.revisions],
        'period': period,
        'components': geometry_data(geometry),
        'measurements': {
            'areas': periodic_area_items(geometry),
            'total_area': periodic_total_area(geometry),
            'boundary': periodic_boundary_measurements(geometry),
        },
        'receipt': {'report_id': report_id, 'batch_id': binding.batch_id,
                    'status': 'published'},
    }
    reconcile_periodic_report(report)
    return report


def build_periodic_report(report_id, batch_id, period):
    from .reading import read_report
    identifier(report_id)
    period = normalize_period(period)
    binding = read_sealed_batch(batch_id)
    existing = publication_kind(report_id)
    if existing is not None:
        stored_format, stored_batch, stored_period = existing
        if stored_format != 2 or stored_batch != batch_id or stored_period != period:
            raise ValueError('report identity is already bound')
        return read_report(report_id)
    current_archive().checkpoint('report.before_compute')
    geometry = build_periodic_geometry(binding.revisions, period)
    report = assemble_periodic_report(report_id, binding, period, geometry)
    publish_report(report)
    return read_report(report_id)
