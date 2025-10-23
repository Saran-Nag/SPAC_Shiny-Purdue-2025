#!/usr/bin/env python3
"""
Script to convert environment.yml to requirements.txt
Handles both conda and pip dependencies
"""

import yaml
import sys
from pathlib import Path


def convert_env_to_requirements(env_file='environment.yml', req_file='requirements.txt'):
    """
    Convert environment.yml to requirements.txt
    
    Args:
        env_file (str): Path to environment.yml file
        req_file (str): Path to output requirements.txt file
    """
    
    # Read environment.yml
    try:
        with open(env_file, 'r') as f:
            env_data = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: {env_file} not found!")
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error parsing {env_file}: {e}")
        sys.exit(1)
    
    requirements = []
    
    # Process conda dependencies
    if 'dependencies' in env_data:
        for dep in env_data['dependencies']:
            if isinstance(dep, str):
                # Skip python version and pip itself
                if dep.startswith('python=') or dep == 'pip':
                    continue
                
                # Convert conda package format to pip format
                # conda uses = for version, pip uses ==
                if '=' in dep and not dep.startswith(('git+', 'http')):
                    package, version = dep.split('=', 1)
                    requirements.append(f"{package}=={version}")
                else:
                    requirements.append(dep)
            
            elif isinstance(dep, dict) and 'pip' in dep:
                # Handle pip dependencies
                pip_deps = dep['pip']
                for pip_dep in pip_deps:
                    requirements.append(pip_dep)
    
    # Remove duplicates while preserving order
    seen = set()
    unique_requirements = []
    for req in requirements:
        # Extract package name for duplicate checking
        pkg_name = req.split('==')[0].split('@')[0].split(' ')[0]
        if pkg_name not in seen:
            seen.add(pkg_name)
            unique_requirements.append(req)
    
    # Write requirements.txt
    try:
        with open(req_file, 'w') as f:
            f.write('\n'.join(unique_requirements))
            f.write('\n')  # Add final newline
        
        print(f"Successfully converted {env_file} to {req_file}")
        print(f"Found {len(unique_requirements)} unique packages")
        
        # Display what was converted
        print("\nPackages included:")
        for req in unique_requirements:
            print(f"  - {req}")
            
    except IOError as e:
        print(f"Error writing {req_file}: {e}")
        sys.exit(1)


def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert environment.yml to requirements.txt')
    parser.add_argument('--env', '-e', default='environment.yml', 
                       help='Input environment.yml file (default: environment.yml)')
    parser.add_argument('--req', '-r', default='requirements.txt',
                       help='Output requirements.txt file (default: requirements.txt)')
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not Path(args.env).exists():
        print(f"Error: Input file {args.env} does not exist!")
        sys.exit(1)
    
    convert_env_to_requirements(args.env, args.req)


if __name__ == '__main__':
    main()