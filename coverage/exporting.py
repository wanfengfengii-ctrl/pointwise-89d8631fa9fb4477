import json
from .reading import read_report
from .result_codec import encode_report

PACKAGE_VERSION_V1 = 1
PACKAGE_VERSION_V2 = 2


def make_package(report):
    version = PACKAGE_VERSION_V2 if report.get('format') == 2 else PACKAGE_VERSION_V1
    return json.dumps({'package_version': version, 'report': json.loads(encode_report(report))},
                      sort_keys=True, separators=(',', ':'), ensure_ascii=False)

def export_report(report_id):
    return make_package(read_report(report_id))
