clear all; clc; close all;

test_result = runtests('tests/test_main.m');
disp(test_result.table)