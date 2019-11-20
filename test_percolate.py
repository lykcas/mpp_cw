import pytest
import os
import subprocess
from test import support


def test_file():
  cmd = 'diff -q ./map8_1.pgm ./map24_oracle.pgm'
  result = subprocess.call(cmd, shell=True)
  assert 0 == result