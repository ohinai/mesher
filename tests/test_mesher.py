#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_mesher
----------------------------------

Tests for `mesher` module.
"""

import unittest

from mesher import mesher


class TestMesher(unittest.TestCase):

    def setUp(self):
        pass

    def test_run_files(self):
        test_files = ["march_2015_1.cmd",
                      "march_2015_2.cmd",
                      "march_2015_3.cmd",]
        for file_name in test_files:
            current_mesh = mesher.Mesher()
            current_mesh.do_load_commands(file_name)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
