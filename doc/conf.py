# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sys
from pathlib import Path

sys.path.insert(0, str(Path("..").resolve()))


# -- Project information -----------------------------------------------------

# General information about the project.
project = "pyphot"
copyright = "2016, M. Fouesneau"
author = "M. Fouesneau"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx_automodapi.automodapi",
    "sphinx.ext.autosectionlabel",
    "sphinx_tabs.tabs",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.todo",
    "matplotlib.sphinxext.plot_directive",
    "sphinx_copybutton",
    "sphinx_mdinclude",
    "sphinx_design",
    "myst_nb",
]
numpydoc_show_class_members = False

# myst_nb configurations
source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}
nb_execution_mode = "off"
myst_enable_extensions = ["dollarmath"]
# auto-generate heading anchors down to this level
myst_heading_anchors = 3

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "**/test_**"]

# hide input prompts in notebooks
nbsphinx_prompt_width = "0"

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_theme_options = {
    "path_to_docs": "doc",
    "repository_url": "https://github.com/mfouesneau/pyphot",
    "repository_branch": "master",
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "colab_url": "https://colab.research.google.com/",
        "notebook_interface": "jupyterlab",
    },
    "use_edit_page_button": False,
    "use_issues_button": True,
    "use_repository_button": True,
    "use_download_button": True,
    "use_source_button": True,
}
html_baseurl = "https://mfouesneau.github.io/pyphot/"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# hide coppy button on outputs
copybutton_selector = "div:not(.output) > div.highlight pre"

# -- Options for autodoc ----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#configuration

# Automatically extract typehints when specified and place them in
# descriptions of the relevant function/method.
autodoc_typehints = "description"

# whether the types of undocumented parameters and return values are documented
autodoc_typehints_description_target = "all"

# controls the format of typehints Suppress the leading module names of the typehints
autodoc_typehints_format = "short"

# Don't show class signature with the class' name.
autodoc_class_signature = "separated"

# remove module name prefix from functions
add_module_names = False

# Eliminate duplicate object warnings while preserving full API documentation coverage.
# https://github.com/dougborg/katana-openapi-client/pull/10
autoapi_python_class_content = "init"

# separate class and init docstrings
autoclass_content = "both"

# Define the order in which automodule and autoclass members are listed
autodoc_member_order = "alphabetical"

autodoc_default_options = {"special-members": "__init__"}

# Prevent imported members from being fully documented at module level
# They will appear as links/references instead
# ! That does not work!
autodoc_inherit_docstrings = False


def autodoc_skip_member(app, what, name, obj, skip, options):
    """
    Skip imported members from submodules to prevent documentation duplication.
    """
    if skip:
        return True

    # Only apply to module-level members
    if what != "module":
        return False

    # Get the documenter from the call stack
    import inspect

    frame = inspect.currentframe()
    try:
        # Walk up the call stack to find the ModuleDocumenter
        while frame:
            local_vars = frame.f_locals
            if "self" in local_vars:
                documenter = local_vars["self"]
                # Check if this is a ModuleDocumenter with the info we need
                if hasattr(documenter, "modname") and hasattr(documenter, "object"):
                    current_module = documenter.modname

                    # Check if obj comes from a submodule
                    if hasattr(obj, "__module__"):
                        obj_module = obj.__module__

                        # Skip if object is from a submodule
                        if obj_module != current_module and obj_module.startswith(
                            current_module + "."
                        ):
                            return True
                    break
            frame = frame.f_back
    finally:
        del frame

    return False


def setup(app):
    app.connect("autodoc-skip-member", autodoc_skip_member)
