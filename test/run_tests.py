#!/usr/bin/env python3
# coding: utf8

import gzip
import hashlib
import itertools
import logging
import os
import tempfile
import subprocess
import sys

logging.basicConfig(level=logging.DEBUG, format='%(message)s')

'''
Assume arboresence like :
|-- HmnTrimmer
|   |-- HmnTrimmer
|   |-- src
|   |-- lib
|   |-- test
|       |-- run_tests.py
|       |   |-- GoldInput
|       |   |-- GoldOutput

'''


def md5ForFile(f, block_size=2**20):
    """Compute MD5 of a file.

    Taken from http://stackoverflow.com/a/1131255/84349.
    """
    md5 = hashlib.md5()
    while True:
        data = f.read(block_size)
        if not data:
            break
        md5.update(data)
    return md5.hexdigest()

class TestConf(object):
    def __init__(self, program, category, name, args, to_diff):
        self.program = program
        self.category = category
        self.name = name
        self.args = args
        self.to_diff = to_diff

    def __repr__(self):
        fmt = '\tTest : %s - %s\n\tArgs : %s\n\tDiff : %s\n'
        return fmt % (self.category, self.name, ' '.join(self.args), ','.join([' '.join(x) for x in self.to_diff]))

    def commandLineArgs(self):
        args = [x for x in self.args if x != '']
        args = [self.program] + args
        return args

def runTest(test_conf):
    def print_error(args, retcode, stdout, stderr):
        logging.error('Return code is %d', retcode)
        if stdout:
            logging.error('--- stdout begin --')
            logging.error(stdout)
            logging.error('--- stdout end --')
        if stderr:
            logging.error('-- stderr begin --')
            logging.error(stderr)
            logging.error('-- stderr end --')
        return False
    # Execute the program.
    #logging.debug('Executing "%s"', ' '.join(test_conf.commandLineArgs()))
    try:
        process = subprocess.Popen(test_conf.commandLineArgs(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if stdout:
            stdout = stdout.decode("utf8")
        if stderr:
            stdeff = stderr.decode("utf8")
        retcode = process.returncode
        if retcode != 0:
            return print_error(test_conf.commandLineArgs(), retcode, stdout, stderr)
    except Exception as e:
        # Print traceback.
        import traceback
        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_exception(exc_type, exc_value, exc_traceback)
        fmt = 'ERROR (when executing "%s"): %s'
        logging.error(fmt % (' '.join(test_conf.commandLineArgs()), e))
        return False

    # Handle error of program, indicated by return code != 0.
    if retcode != 0:
        return print_error(test_conf.commandLineArgs(), retcode, stdout, stderr)

    # Compare results with expected results, if the expected and actual result
    # are not equal then print diffs.
    results = []
    for tuple_ in test_conf.to_diff:
        expected_path, result_path = tuple_[:2]
        mode = 'binary'
        if len(tuple_) >= 3:
            if tuple_[2] == 'md5':
                mode = 'binary'
            elif tuple_[2] == 'gzip':
                mode = 'gzip'
            elif tuple_[2] == 'empty':
                mode = 'empty'
        try:
            if mode == 'gzip':
                f = gzip.open(expected_path, 'rb')
                expected_md5 = md5ForFile(f)
                f.close()
                f = gzip.open(result_path, 'rb')
                result_md5 = md5ForFile(f)
                f.close()
                if expected_md5 == result_md5:
                    results.append(0)
                else:
                    tpl = (expected_path, expected_md5, result_md5, result_path)
                    logging.error('md5(gzip(%s)) == %s != %s == md5(gzip(%s))' % tpl)
                    return False
            elif mode == 'binary':
                with open(expected_path, 'rb') as f:
                    expected_md5 = md5ForFile(f)
                with open(result_path, 'rb') as f:
                    result_md5 = md5ForFile(f)
                if expected_md5 == result_md5:
                    results.append(0)
                else:
                    tpl = (expected_path, expected_md5, result_md5, result_path)
                    logging.error('md5(%s) == %s != %s == md5(%s)' % tpl)
                    return False
            elif mode == 'empty':
                if os.stat(result_path).st_size == 0:
                    results.append(0)
                else:
                    logging.error('File %s is not empty'%(result_path,))
                    return False        
        except Exception as e:
            fmt = 'Error when trying to compare %s to %s: %s ' + str(type(e))
            logging.error(fmt % (expected_path, result_path, e))
            return False
    if sum(results) > 0:
        return False
    return True

def clean_up_files(files):
    for filename in files:
        os.remove(filename)
    files = []

def create_tmp_files(temp_files, files, nb=1, exts=None):
    files.clear()
    if type(exts) == str:
        exts = [exts] * nb
    for x, ext in zip(range(nb), exts):
        if ext:
            temp = tempfile.mkstemp(suffix=ext)
        else:
            temp = tempfile.mkstemp()

        files.append(temp[1])
    temp_files += files

def main():
    # Defined paths.
    dir_test = os.path.dirname(os.path.realpath(__file__))
    path_program = os.path.join(os.path.dirname(dir_test), 'HmnTrimmer')
    path_gold_input = os.path.join(dir_test, 'GoldInput')
    path_gold_output = os.path.join(dir_test, 'GoldOutput')

    # Init var.
    conf_list = []
    TMPFILES, temp_files = [], [] # First is global, second is purposed for each test

    # ============================================================
    # GenTrim.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 1, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenTrim',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--length-min', '10'],
        to_diff=[(os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
              temp_files[0])]
    )
    conf_list.append(conf)

    # B.
    create_tmp_files(TMPFILES, temp_files, 1, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenTrim',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--length-min', '21'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-A.R1.fastq'), temp_files[0])]
    )
    conf_list.append(conf)

    # C.
    create_tmp_files(TMPFILES, temp_files, 1, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenTrim',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--length-min', '50'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-B.R1.single.fastq'), temp_files[0])]
    )
    conf_list.append(conf)

    # D.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenTrim',
        name = 'D',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '10'],
        to_diff=[(os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # E.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenTrim',
        name = 'E',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '20'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-A.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'LENGTHMIN-A.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # F.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenTrim',
        name = 'F',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '50'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-B.R1.paired.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'LENGTHMIN-B.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # ============================================================
    # GenDiscard.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenDiscard',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--output-fastq-discard', temp_files[0],
            '--output-fastq-forward', temp_files[1],
            '--length-min', '10'],
        to_diff=[('', temp_files[0], 'empty')]
    )
    conf_list.append(conf)

    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenDiscard',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--output-fastq-discard', temp_files[0],
            '--output-fastq-forward', temp_files[1],
            '--length-min', '20'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-A.R1.single.discard.fastq'), temp_files[0])]
    )
    conf_list.append(conf)

    # C.
    create_tmp_files(TMPFILES, temp_files, 3, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenDiscard',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-discard', temp_files[0],
            '--output-fastq-forward', temp_files[1],
            '--output-fastq-reverse', temp_files[2],
            '--length-min', '50'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-B.discard.fastq'), temp_files[0])]
    )
    conf_list.append(conf)

    # ============================================================
    # GenFormat.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 1, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenFormat',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-interleaved', temp_files[0],
            '--length-min', '10'],
        to_diff=[(os.path.join(path_gold_input, 'LENGTHMIN.Interleaved.fastq'), temp_files[0])]
    )
    conf_list.append(conf)

    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenFormat',
        name = 'B',
        args = ['--input-fastq-interleaved', os.path.join(path_gold_input, 'LENGTHMIN.Interleaved.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '10'],
        to_diff=[(os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # C.
    create_tmp_files(TMPFILES, temp_files, 1, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenFormat',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-interleaved', temp_files[0],
            '--length-min', '50'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-B.Interleaved.fastq'), temp_files[0])]
    )
    conf_list.append(conf)

    # D.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenFormat',
        name = 'D',
        args = ['--input-fastq-interleaved', os.path.join(path_gold_input, 'LENGTHMIN.Interleaved.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '50'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-B.R1.paired.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'LENGTHMIN-B.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # ============================================================
    # GenThread.
    # ============================================================
    tests = [ (chr(name), str(thread), '10') for name, thread in zip(range(ord('A'), ord('E')), range(2,10,2)) ]
    for test in tests:
        create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
        conf = TestConf(
            program = path_program,
            category = 'GenThread',
            name = test[0],
            args = ['--input-fastq-forward', os.path.join(path_gold_input, 'BIG.R1.fastq'),
                '--input-fastq-reverse', os.path.join(path_gold_input, 'BIG.R2.fastq'),
                '--output-fastq-forward', temp_files[0],
                '--output-fastq-reverse', temp_files[1],
                '--reads-batch', '100',
                '--threads', test[1],
                '--length-min', test[2]],
            to_diff=[(os.path.join(path_gold_input, 'BIG.R1.fastq'), temp_files[0]),
                (os.path.join(path_gold_input, 'BIG.R2.fastq'), temp_files[1])]
        )
        conf_list.append(conf)

    tests = [ (chr(name), str(thread), '55') for name, thread in zip(range(ord('E'), ord('I')), range(2,10,2)) ]
    for test in tests:
        create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
        conf = TestConf(
            program = path_program,
            category = 'GenThread',
            name = test[0],
            args = ['--input-fastq-forward', os.path.join(path_gold_input, 'BIG.R1.fastq'),
                '--input-fastq-reverse', os.path.join(path_gold_input, 'BIG.R2.fastq'),
                '--output-fastq-forward', temp_files[0],
                '--output-fastq-reverse', temp_files[1],
                '--reads-batch', '100',
                '--threads', test[1],
                '--length-min', test[2]],
            to_diff=[(os.path.join(path_gold_output, 'BIG-B.R1.fastq'), temp_files[0]),
                (os.path.join(path_gold_output, 'BIG-B.R2.fastq'), temp_files[1])]
        )
        conf_list.append(conf)
    # ============================================================
    # GenCompress.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq.gz")
    conf = TestConf(
        program = path_program,
        category = 'GenCompress',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'BIG.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'BIG.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '10'],
        to_diff=[(os.path.join(path_gold_input, 'BIG.R1.fastq.gz'), temp_files[0], 'gzip'),
            (os.path.join(path_gold_input, 'BIG.R2.fastq.gz'), temp_files[1], 'gzip')]
    )
    conf_list.append(conf)
    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'GenCompress',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'BIG.R1.fastq.gz'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'BIG.R2.fastq.gz'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '10'],
        to_diff=[(os.path.join(path_gold_input, 'BIG.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_input, 'BIG.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # C.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq.gz")
    conf = TestConf(
        program = path_program,
        category = 'GenCompress',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'BIG.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'BIG.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '55'],
        to_diff=[(os.path.join(path_gold_output, 'BIG-B.R1.fastq.gz'), temp_files[0], 'gzip'),
            (os.path.join(path_gold_output, 'BIG-B.R2.fastq.gz'), temp_files[1], 'gzip')]
    )
    conf_list.append(conf)

    # ============================================================
    # TrimLengthMin.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimLengthMin',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '10'],
        to_diff=[(os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimLengthMin',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '50'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-B.R1.paired.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'LENGTHMIN-B.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # C.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimLengthMin',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'LENGTHMIN.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'LENGTHMIN.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--length-min', '100'],
        to_diff=[(os.path.join(path_gold_output, 'LENGTHMIN-C.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'LENGTHMIN-C.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # ============================================================
    # TrimQualSld.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualSld',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALSLD.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALSLD.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-sliding-window', '1:5'],
        to_diff=[(os.path.join(path_gold_input, 'QUALSLD.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_input, 'QUALSLD.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualSld',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALSLD.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALSLD.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-sliding-window', '20:25'],
        to_diff=[(os.path.join(path_gold_output, 'QUALSLD-A.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALSLD-A.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # C.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualSld',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALSLD.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALSLD.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-sliding-window', '10:10'],
        to_diff=[(os.path.join(path_gold_output, 'QUALSLD-B.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALSLD-B.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # D.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualSld',
        name = 'D',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALSLD.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALSLD.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-sliding-window', '17:4'],
        to_diff=[(os.path.join(path_gold_output, 'QUALSLD-C.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALSLD-C.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # E.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualSld',
        name = 'E',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALSLD.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALSLD.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-sliding-window', '17:10'],
        to_diff=[(os.path.join(path_gold_output, 'QUALSLD-D.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALSLD-D.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # F.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualSld',
        name = 'F',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALSLD.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALSLD.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-sliding-window', '4:1'],
        to_diff=[(os.path.join(path_gold_output, 'QUALSLD-E.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALSLD-E.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # G.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualSld',
        name = 'G',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALSLD.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALSLD.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-sliding-window', '4:4'],
        to_diff=[(os.path.join(path_gold_output, 'QUALSLD-F.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALSLD-F.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # ============================================================
    # TrimQualTail.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualTail',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALTAIL.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALTAIL.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-tail', '2:2'],
        to_diff=[(os.path.join(path_gold_output, 'QUALTAIL-A.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALTAIL-A.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualTail',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALTAIL.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALTAIL.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-tail', '5:2'],
        to_diff=[(os.path.join(path_gold_output, 'QUALTAIL-B.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALTAIL-B.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # C.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualTail',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALTAIL.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALTAIL.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-tail', '5:2:60'],
        to_diff=[(os.path.join(path_gold_output, 'QUALTAIL-C.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALTAIL-C.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # D.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualTail',
        name = 'D',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALTAIL.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALTAIL.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-tail', '5:5:60'],
        to_diff=[(os.path.join(path_gold_output, 'QUALTAIL-C.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALTAIL-C.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # E.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualTail',
        name = 'E',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALTAIL.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALTAIL.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-tail', '5:5:70'],
        to_diff=[(os.path.join(path_gold_output, 'QUALTAIL-D.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALTAIL-D.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # F.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualTail',
        name = 'F',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALTAIL.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALTAIL.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-tail', '3:10:46'],
        to_diff=[(os.path.join(path_gold_output, 'QUALTAIL-E.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALTAIL-E.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # G.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimQualTail',
        name = 'G',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'QUALTAIL.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'QUALTAIL.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--quality-tail', '3:10:81'],
        to_diff=[(os.path.join(path_gold_output, 'QUALTAIL-F.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'QUALTAIL-F.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # ============================================================
    # TrimInfoDust.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimInfoDust',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'INFODUST.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'INFODUST.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--information-dust', '7'],
        to_diff=[(os.path.join(path_gold_input, 'INFODUST.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_input, 'INFODUST.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimInfoDust',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'INFODUST.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'INFODUST.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--information-dust', '3'],
        to_diff=[(os.path.join(path_gold_output, 'INFODUST-A.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'INFODUST-A.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)
    # C.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimInfoDust',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'INFODUST.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'INFODUST.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--information-dust', '5'],
        to_diff=[(os.path.join(path_gold_output, 'INFODUST-B.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'INFODUST-B.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # ============================================================
    # TrimInfoN.
    # ============================================================
    # A.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimN',
        name = 'A',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'INFON.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'INFON.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--information-n', '5'],
        to_diff=[(os.path.join(path_gold_input, 'INFON.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_input, 'INFON.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # B.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimN',
        name = 'B',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'INFON.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'INFON.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--information-n', '1'],
        to_diff=[(os.path.join(path_gold_output, 'INFON-A.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'INFON-A.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # C.
    create_tmp_files(TMPFILES, temp_files, 2, ".fastq")
    conf = TestConf(
        program = path_program,
        category = 'TrimN',
        name = 'C',
        args = ['--input-fastq-forward', os.path.join(path_gold_input, 'INFON.R1.fastq'),
            '--input-fastq-reverse', os.path.join(path_gold_input, 'INFON.R2.fastq'),
            '--output-fastq-forward', temp_files[0],
            '--output-fastq-reverse', temp_files[1],
            '--information-n', '4'],
        to_diff=[(os.path.join(path_gold_output, 'INFON-B.R1.fastq'), temp_files[0]),
            (os.path.join(path_gold_output, 'INFON-B.R2.fastq'), temp_files[1])]
    )
    conf_list.append(conf)

    # ============================================================
    # Execute the tests.
    # ============================================================
    failures = []
    for conf in conf_list:
        res = runTest(conf)
        if res:
            logging.info('Test %s %s\tOK'%(conf.category, conf.name,))
        else:
            failures.append(conf)

    if len(failures) > 0:
        logging.info('==============================')
        for conf in failures:
            logging.error('Test %s %s\tFAILED'%(conf.category, conf.name,))
            logging.debug('RunTest :\n%s\n%s'%(' '.join([conf.program] + conf.args), conf))    
        logging.info('==============================')
        
    # Cleanup.
    logging.debug('Clean up.')
    clean_up_files(TMPFILES)

    logging.info('==============================')
    logging.info('     total tests: %s' % (len(conf_list),))
    logging.info('    failed tests: %s' % (len(failures),))
    logging.info('successful tests: %s' % (len(conf_list) - len(failures),))
    logging.info('==============================')
    # Compute and return return code.
    return len(failures) != 0

if __name__ == '__main__':
    sys.exit(main())
