from .result_codec import parse_json, unpack, validate_report_structure
from .coverage_union import require_keys
from .reconciliation import reconcile_report

SUPPORTED_PACKAGE_VERSIONS = (1, 2)

def parse_package(payload):
    raw = parse_json(payload)
    require_keys(raw, ('package_version', 'report'))
    if type(raw['package_version']) is not int or isinstance(raw['package_version'], bool) \
            or raw['package_version'] not in SUPPORTED_PACKAGE_VERSIONS:
        raise ValueError('unsupported package version')
    report = unpack(raw['report'])
    validate_report_structure(report)
    if report['format'] != raw['package_version']:
        raise ValueError('package version does not match report format')
    return report

def open_package(payload):
    report = parse_package(payload)
    reconcile_report(report)
    return report
