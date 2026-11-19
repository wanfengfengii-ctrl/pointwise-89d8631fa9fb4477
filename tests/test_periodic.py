import json
import sqlite3
import tempfile
import unittest
from fractions import Fraction as F
from pathlib import Path

import coverage as api
from coverage import open_archive, open_package, put_region, seal_batch

P = F(10)


def box(a, b, c, d):
    return [(a, b), (c, b), (c, d), (a, d)]


def seal(shapes, batch='day', shift=0):
    members = []
    for name, outer, holes in shapes:
        def moved(ring):
            return [(F(x) + shift, F(y)) for x, y in ring]
        revision = put_region(name, moved(outer), [moved(h) for h in holes])
        members.append((name, revision))
    seal_batch(batch, members)
    return batch


class TestPeriodicIntake(unittest.TestCase):
    def test_period_validation(self):
        with open_archive(':memory:'):
            seal([('A', box(9, 0, 11, 1), [])])
            for invalid in (0, -3, True, 10.0, '10', None, F(0)):
                with self.subTest(invalid=invalid):
                    with self.assertRaises(ValueError):
                        api.build_periodic_report('r', 'day', invalid)

    def test_span_and_edge_limits(self):
        with open_archive(':memory:'):
            put_region('wide', box(0, 0, 10, 1))
            seal_batch('span', [('wide', 1)])
            with self.assertRaises(ValueError):
                api.build_periodic_report('r', 'span', P)
            many = [('R%d' % i, box(i, 0, i + 1, 1), []) for i in range(21)]
            seal(many, 'many')
            with self.assertRaises(ValueError):
                api.build_periodic_report('big', 'many', P)

    def test_missing_batch_is_key_error(self):
        with open_archive(':memory:'):
            with self.assertRaises(KeyError):
                api.build_periodic_report('r', 'absent', P)


class TestPeriodicReports(unittest.TestCase):
    def build(self, shapes, report='r', batch='day', period=P, shift=0):
        seal(shapes, batch, shift)
        return api.build_periodic_report(report, batch, period)

    def test_shape_and_freeze_rules(self):
        with open_archive(':memory:'):
            report = self.build([('A', box(9, 0, 11, 1), [])])
            self.assertEqual(set(report), {'format', 'report_id', 'batch_id', 'binding',
                                           'inputs', 'period', 'components', 'measurements',
                                           'receipt'})
            self.assertEqual(report['format'], 2)
            self.assertEqual(report['period'], P)
            self.assertIs(type(report['period']), F)
            self.assertEqual(api.build_periodic_report('r', 'day', 10), report)
            with self.assertRaises(ValueError):
                api.build_periodic_report('r', 'day', F(11))
            api.build_report('flat', 'day')
            with self.assertRaises(ValueError):
                api.build_report('r', 'day')
            with self.assertRaises(ValueError):
                api.build_periodic_report('flat', 'day', P)
            payload = json.loads(api.export_report('r'))
            self.assertEqual(payload['package_version'], 2)
            self.assertEqual(payload['report']['format'], 2)

    def test_winding_band(self):
        with open_archive(':memory:'):
            report = self.build([('A', box(0, 0, 6, 1), []),
                                 ('B', box(5, 0, 11, 1), [])], report='band')
            self.assertEqual(report['measurements']['total_area'], F(10))
            windings = sorted(b['winding'] for c in report['components'] for b in c['boundaries'])
            self.assertEqual(windings, [-1, 1])
            for component in report['components']:
                for boundary in component['boundaries']:
                    for segment in boundary['segments']:
                        self.assertEqual(segment['delta'][1], 0)
                        self.assertTrue(segment['sources'])

    def test_point_versus_edge_contact(self):
        with open_archive(':memory:'):
            seal([('A', box(9, 0, 10, 1), []),
                  ('B', box(0, 1, 1, 2), [])], 'point')
            point = api.build_periodic_report('point', 'point', P)
            self.assertEqual(len(point['components']), 2)
            self.assertEqual(point['measurements']['total_area'], F(2))
        with open_archive(':memory:'):
            seal([('A', box(9, 0, 10, 1), []),
                  ('B', box(0, 0, 1, 1), [])], 'edge')
            edge = api.build_periodic_report('edge', 'edge', P)
            self.assertEqual(len(edge['components']), 1)
            self.assertEqual(edge['measurements']['total_area'], F(2))

    def test_cross_seam_hole(self):
        with open_archive(':memory:'):
            report = self.build([('A', box(7, 0, 15, 4), [box(9, 1, 11, 3)]),
                                 ('B', box(4, 0, 8, 4), [])], report='hole')
            self.assertEqual(report['measurements']['total_area'], F(36))
            windings = sorted(b['winding'] for c in report['components'] for b in c['boundaries'])
            self.assertEqual(windings, [-1, 0, 1])
            hole_area = F(4)
            self.assertEqual(report['measurements']['areas'][0]['area'], F(36))
            self.assertTrue(all(item['area'] > 0 for item in report['measurements']['areas']))

    def test_duplicate_coincident_sources(self):
        with open_archive(':memory:'):
            report = self.build([('A', box(9, 0, 11, 1), []),
                                 ('B', list(reversed(box(-1, 0, 1, 1))), [])],
                                report='dup')
            for component in report['components']:
                for boundary in component['boundaries']:
                    for segment in boundary['segments']:
                        self.assertEqual({s['region_id'] for s in segment['sources']}, {'A', 'B'})

    def test_translation_equivalence(self):
        shapes = [('A', box(9, 0, 11, 1), [])]
        with open_archive(':memory:'):
            first = self.build(shapes, report='r', shift=0)
        with open_archive(':memory:'):
            second = self.build(shapes, report='r', shift=-20)
        self.assertEqual(first['measurements'], second['measurements'])
        first_windings = [b['winding'] for c in first['components'] for b in c['boundaries']]
        second_windings = [b['winding'] for c in second['components'] for b in c['boundaries']]
        self.assertEqual(first_windings, second_windings)

    def test_fractional_slanted_edge(self):
        with open_archive(':memory:'):
            report = self.build([('T', [(9, 0), (11, 0), (10, 2)], [])], report='slant')
            self.assertEqual(report['measurements']['total_area'], F(2))
            boundary = report['components'][0]['boundaries'][0]
            bid = boundary['boundary_id']
            sample = api.sample_boundary('slant', bid, 0, F(1, 2))
            self.assertEqual(sample['point'][1], F(1))

    def test_length_squared_terms(self):
        with open_archive(':memory:'):
            report = self.build([('A', box(9, 0, 11, 1), [])])
            walk = report['components'][0]['boundaries'][0]
            expected = [s['delta'][0] ** 2 + s['delta'][1] ** 2 for s in walk['segments']]
            measured = next(item['length_squared']
                            for item in report['measurements']['boundary']
                            if item['boundary_id'] == walk['boundary_id'])
            self.assertEqual(measured, expected)


class TestPeriodicConsumption(unittest.TestCase):
    def test_reopen_and_offline_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'periodic.sqlite'
            with open_archive(path):
                seal([('A', box(0, 0, 6, 1), []), ('B', box(5, 0, 11, 1), [])])
                report = api.build_periodic_report('band', 'day', P)
                saved = api.export_report('band')
                bid = report['components'][0]['boundaries'][0]['boundary_id']
                sample = api.sample_boundary('band', bid, 0, F(1, 3))
            with open_archive(path):
                self.assertEqual(api.export_report('band'), saved)
                self.assertEqual(api.sample_boundary('band', bid, 0, F(1, 3)), sample)
                self.assertEqual(api.read_report('band'), report)
            self.assertEqual(open_package(saved), report)

    def test_package_rejects_damage(self):
        with open_archive(':memory:'):
            seal([('A', box(9, 0, 11, 1), [])])
            api.build_periodic_report('r', 'day', P)
            payload = json.loads(api.export_report('r'))
            for damage in (
                lambda p: p.update({'package_version': 1}),
                lambda p: p['report'].update({'format': 1}),
                lambda p: p['report']['components'][0]['boundaries'][0].update({'winding': 7}),
                lambda p: p['report']['components'][0]['boundaries'][0]['segments'][0]['sources'].pop(),
                lambda p: p['report']['measurements']['areas'][0].update({'area': {'q': [99, 1]}}),
            ):
                broken = json.loads(json.dumps(payload))
                damage(broken)
                with self.subTest(damage=damage):
                    with self.assertRaises(ValueError):
                        open_package(json.dumps(broken))

    def test_sample_parameter_and_identity_errors(self):
        with open_archive(':memory:'):
            seal([('A', box(9, 0, 11, 1), [])])
            report = api.build_periodic_report('r', 'day', P)
            bid = report['components'][0]['boundaries'][0]['boundary_id']
            count = len(report['components'][0]['boundaries'][0]['segments'])
            for bad in (0, F(-1), F(2), True, 0.5):
                with self.assertRaises(ValueError):
                    api.sample_boundary('r', bid, 0, bad)
            with self.assertRaises(KeyError):
                api.sample_boundary('absent', bid, 0, F(0))
            with self.assertRaises(KeyError):
                api.sample_boundary('r', 'missing', 0, F(0))
            with self.assertRaises(KeyError):
                api.sample_boundary('r', bid, count, F(0))
            point = api.sample_boundary('r', bid, count - 1, F(1))
            self.assertTrue(0 <= point['point'][0] < P)

    def test_publication_failure_leaves_nothing(self):
        with open_archive(':memory:') as archive:
            seal([('A', box(9, 0, 11, 1), [])])
            def fail(name, connection):
                if name == 'report.after_payload':
                    connection.execute('INSERT INTO nope VALUES(1)')
            archive.failure_hook = fail
            with self.assertRaises(sqlite3.OperationalError):
                api.build_periodic_report('r', 'day', P)
            for table in ('reports', 'receipts'):
                self.assertEqual(archive.connection.execute(
                    'SELECT count(*) FROM ' + table + " WHERE report_id='r'").fetchone()[0], 0)
            archive.failure_hook = None
            self.assertEqual(api.build_periodic_report('r', 'day', P)
                             ['measurements']['total_area'], F(2))

    def test_upstream_correction_does_not_drift(self):
        with open_archive(':memory:'):
            seal([('A', box(9, 0, 11, 1), [])])
            report = api.build_periodic_report('r', 'day', P)
            saved = api.export_report('r')
            report['components'].clear()
            put_region('A', box(0, 0, 5, 5))
            self.assertEqual(api.export_report('r'), saved)
            self.assertEqual(api.read_report('r')['measurements']['total_area'], F(2))


if __name__ == '__main__':
    unittest.main()
