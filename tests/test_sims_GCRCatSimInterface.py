"""
Example unit tests for sims_GCRCatSimInterface package
"""
import unittest
import desc.sims_gcrcatsiminterface

class sims_GCRCatSimInterfaceTestCase(unittest.TestCase):
    def setUp(self):
        self.message = 'Hello, world'

    def tearDown(self):
        pass

    def test_run(self):
        foo = desc.sims_gcrcatsiminterface.sims_GCRCatSimInterface(self.message)
        self.assertEquals(foo.run(), self.message)

    def test_failure(self):
        self.assertRaises(TypeError, desc.sims_gcrcatsiminterface.sims_GCRCatSimInterface)
        foo = desc.sims_gcrcatsiminterface.sims_GCRCatSimInterface(self.message)
        self.assertRaises(RuntimeError, foo.run, True)

if __name__ == '__main__':
    unittest.main()
