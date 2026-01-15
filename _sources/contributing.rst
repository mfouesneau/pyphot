Contributing
============

Please open a new issue or new pull request for bugs, feedback, or new features you would like to see. If there is an issue you would like to work on, please leave a comment, and we will be happy to assist. New contributions and contributors are very welcome!

This document is a guideline to help you get started with contributing to PyPhot. It covers the following topics:

* `Setting up your development environment <#setting-up-your-development-environment>`_
* `Running tests <#testing>`_
* `Style guide for writing code <#style-guide>`_
* `Documentation <#contributing-to-the-documentation>`_

Please also have a look at our `code of conduct <https://github.com/mfouesneau/pyphot/blob/master/CODE_OF_CONDUCT.md>`_.

Setting up your development environment
---------------------------------------

Virtual environment
~~~~~~~~~~~~~~~~~~~
We recommend using a virtual environment to isolate your development environment from your system's Python installation. You can create a virtual environment using the following command:

* using Python's built-in venv module

.. code-block:: bash

    python -m venv .venv
    source .venv/bin/activate

* using `uv <https://docs.astral.sh/uv/>`_

.. code-block:: bash

    # install uv
    curl -LsSf https://astral.sh/uv/install.sh | sh
    # create a virtual environment
    uv venv
    # activate the virtual environment
    source .venv/bin/activate

You can then clone the repository and install it in editable mode with the extra development dependencies.

.. code-block:: bash

    git clone https://github.com/mfouesneau/pyphot.git
    cd pyphot
    pip install -e .[ci,testing,docs]

* `ci` contains the linters and pre-commit tools.

* `testing` contains the packages to run the unit tests and generate test coverage reports.

* `docs` contains the package to compile the documentation with sphinx and myst.

You may not want to install all the dependencies depending on your contribution.

Setting up pre-commit hooks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This step is optional but highly recommended if you want to ensure your code will not be stuck with a failing Pull Request. You can iterate on your code locally and avoid unnecessary pushes.

Once you have developed your new functionality, you'll want to commit it to the repo. We employ a pre-commit hooks to ensure any code you commit. These hooks will guard against against committing of large files, sanitise your scripts and Jupyter notebooks, and maybe most importantly, will run ruff in both linter and formatter mode to ensure your code is formatted correctly.

This requires a small amount of set-up on your part, some of which was done when you installed the optional `.[ci]` dependencies above. The rest of the setup requires you run

.. code-block:: bash

    prek install

at the root of the repo to activate the pre-commit hooks. We use `prek <https://prek.j178.dev/>`_ which is a faster alternative to `pre-commit <https://pre-commit.com/>`_.

If you would like to test whether it works you can run

.. code-block:: bash

    prek run --all-files

to run the pre-commit hook on the whole repo. You should see each stage complete without issue in a clean state.

Testing
-------
Unit tests are a crucial part of the development process. They help ensure that your code works as expected and that any changes you make don't break existing functionality.
These are not included in the pre-commit hooks (as they take a substantial amount of time to run), but you can run them manually.

To run the unit tests, you can use the following command:

.. code-block:: bash

    pytest

This will execute all the tests in the `unittests` directory.

GitHub workflow will automatically run linters (Ruff) and tests on any contributions; builds that fail these tests will not be accepted.

Style guide
-----------

All new PRs should follow these guidelines. We adhere to the PEP-8 style guide, and as described above this is verified with ruff.

The ruff configuration is defined in our pyproject.toml so there's no need to configure it yourself, we've made all the decisions for you (for better or worse). Any merge request will be checked with the ruff linter and must pass before being eligible to merge.

We use the `Numpy style guide <https://numpydoc.readthedocs.io/en/latest/format.html>`_ for docstrings.

Contributing to the Documentation
---------------------------------

The documentation is written in a combination of reStructuredText, Markdown, Jupyter notebooks and Python scripts. Adding content should be relatively simple if you follow the instructions below.

We use `Sphinx <https://www.sphinx-doc.org/en/master/>`_ for documentation and `myst <https://myst-parser.readthedocs.io/en/latest/>`_ for Markdown and Notebook rendering. The necessary extensions are installed with `.[docs]`.

The published documentation reflects the current distribution available on PyPI. If you would like to see the current development version in your branch or the main branch, you will have to build the documentation locally. To do so, navigate to the docs directory and run:

The documentation is built using the following command (from the `doc` directory):

.. code-block:: bash

    make clean && make html

This will generate the HTML documentation in the `_build/html` directory.

You can serve your version locally with

.. code-block:: bash

    python -m http.server -d _build/html

This will start a local server (typically at http://localhost:8000).

GitHub workflow will not automatically build the documentation on any contributions; It only runs on the main branch.
