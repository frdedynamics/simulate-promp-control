clear all; clc; close all;

test_result = runtests('IncludeSubfolders', true);
disp(test_result.table)