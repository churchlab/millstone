"""
Convenience module for debugging scripts.
"""

# Since this script is intended to be used from the terminal, setup the
# environment first so that django and model imports work.
from util import setup_django_env
setup_django_env()


def main():
    pass

if __name__ == '__main__':
    main()
