#!/usr/bin/env python
# coding: utf-8

# # Bioinformatics
# 
# This is a jupyterbook compilation of Notebooks for bioinformatics in Python.
# 
# * Biopython

# In[1]:


import IPython
print( IPython.sys_info() )


# ## How to build

# ### Building the book
# 
# If you'd like to develop on and build the bioinformatics book, you should:
# 
# - Clone this repository and run
# - Run `pip install -r requirements.txt` (it is recommended you do this within a virtual environment)
# - (Recommended) Remove the existing `bioinformatics/_build/` directory
# - Run `jupyter-book build bioinformatics/`
# 
# A fully-rendered HTML version of the book will be built in `bioinformatics/_build/html/`.
# 
# ### Hosting the book
# 
# The html version of the book is hosted on the `gh-pages` branch of this repo. A GitHub actions workflow has been created that automatically builds and pushes the book to this branch on a push or pull request to main.
# 
# If you wish to disable this automation, you may remove the GitHub actions workflow and build the book manually by:
# 
# - Navigating to your local build; and running,
# - `ghp-import -n -p -f bioinformatics/_build/html`
# 
# This will automatically push your build to the `gh-pages` branch. More information on this hosting process can be found [here](https://jupyterbook.org/publish/gh-pages.html#manually-host-your-book-with-github-pages).

# In[ ]:





# 
# ```{toctree}
# :hidden:
# :titlesonly:
# :caption: Biopython
# 
# notebooks/biopython/02-Quick-Start.ipynb
# notebooks/biopython/03-Sequence-Objects.ipynb
# notebooks/biopython/04-Sequence-Annotation-Objects.ipynb
# ```
# 
