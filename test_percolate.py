import pytest
import subprocess


def test_file():
  cmd = 'diff -q ./map288.pgm ./map288_oracle.pgm'
  result = subprocess.call(cmd, shell=True)
  assert 0 == result

# def test_check():
#   cmd = 'diff -q ./log ./log_oracle'
#   result = subprocess.call(cmd, shell=True)
#   assert 0 == result