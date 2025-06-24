#!/usr/bin/env python3
"""
Test script for GenPPI Python interface
"""

import os
import sys
from pathlib import Path
from genppi_py.genppi import get_executable_path, run_genppi

def test_executable_path():
    """Test that the executable path is correctly determined"""
    try:
        path = get_executable_path(use_db=False)
        print(f"GenPPI executable path: {path}")
        assert path.exists(), f"Executable not found: {path}"
        
        path_db = get_executable_path(use_db=True)
        print(f"GenPPIDB executable path: {path_db}")
        assert path_db.exists(), f"Executable not found: {path_db}"
        
        print("Executable path test passed!")
        return True
    except Exception as e:
        print(f"Error in executable path test: {e}")
        return False

def test_help_command():
    """Test running the help command"""
    try:
        ret_code = run_genppi(['-help'], use_db=False)
        assert ret_code == 0, f"Help command failed with return code {ret_code}"
        
        ret_code_db = run_genppi(['-help'], use_db=True)
        assert ret_code_db == 0, f"Help command failed with return code {ret_code_db}"
        
        print("Help command test passed!")
        return True
    except Exception as e:
        print(f"Error in help command test: {e}")
        return False

def test_ppmethod_parameter():
    """Test the ppmethod parameter mapping"""
    try:
        from genppi_py.genppi import main
        import sys
        from unittest.mock import patch
        
        # Test ppmethod 1 (should map to -ppcomplete)
        test_args = ['-ppmethod', '1', '-dir', 'test']
        with patch.object(sys, 'argv', ['genppi'] + test_args):
            with patch('genppi_py.genppi.run_genppi') as mock_run:
                mock_run.return_value = 0
                main()
                # Check if -ppcomplete was added to the command args
                call_args = mock_run.call_args[0][0]
                assert '-ppcomplete' in call_args, f"Expected -ppcomplete in args, got: {call_args}"
        
        print("ppmethod parameter test passed!")
        return True
    except Exception as e:
        print(f"Error in ppmethod parameter test: {e}")
        return False

def test_cn_parameter_removed():
    """Test that the -cn parameter was correctly removed"""
    try:
        from genppi_py.genppi import main
        import sys
        from unittest.mock import patch
        
        # Test that -cn parameter is no longer accepted
        test_args = ['-cn', '-dir', 'test']
        with patch.object(sys, 'argv', ['genppi'] + test_args):
            with patch('genppi_py.genppi.run_genppi') as mock_run:
                mock_run.return_value = 0
                try:
                    main()
                    # If we get here, -cn was incorrectly accepted
                    assert False, "-cn parameter should have been removed"
                except SystemExit:
                    # This is expected - argparse should exit with error for unknown argument
                    pass
        
        print("CN parameter removal test passed!")
        return True
    except Exception as e:
        print(f"Error in CN parameter removal test: {e}")
        return False

def test_default_behavior():
    """Test that default behavior works without explicit -pp flag"""
    try:
        from genppi_py.genppi import main
        import sys
        from unittest.mock import patch
        
        # Test basic usage without -pp flag
        test_args = ['-dir', 'test']
        with patch.object(sys, 'argv', ['genppi'] + test_args):
            with patch('genppi_py.genppi.run_genppi') as mock_run:
                mock_run.return_value = 0
                main()
                # Check that the command was called (meaning no errors in argument parsing)
                call_args = mock_run.call_args[0][0]
                assert '-dir' in call_args, f"Expected -dir in args, got: {call_args}"
                # Check that the directory argument is present (may be expanded to absolute path)
                dir_index = call_args.index('-dir')
                assert dir_index + 1 < len(call_args), f"Expected directory argument after -dir"
                dir_arg = call_args[dir_index + 1]
                assert 'test' in dir_arg, f"Expected 'test' in directory argument, got: {dir_arg}"
        
        print("Default behavior test passed!")
        return True
    except Exception as e:
        print(f"Error in default behavior test: {e}")
        return False

if __name__ == '__main__':
    success = True
    success = test_executable_path() and success
    success = test_help_command() and success
    success = test_ppmethod_parameter() and success
    success = test_cn_parameter_removed() and success
    success = test_default_behavior() and success
    
    if success:
        print("All tests passed!")
        sys.exit(0)
    else:
        print("Some tests failed!")
        sys.exit(1)