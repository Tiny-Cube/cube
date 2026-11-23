project = 'Wren'
from pathlib import Path
import jieba  # Required: Sphinx otherwise silently omits Chinese search terms.
author = 'Tiny-Cube'
copyright = '2026, Tiny-Cube'
release = '0.2.0'
extensions = ['myst_parser', 'sphinx.ext.mathjax']
myst_enable_extensions = ['dollarmath', 'colon_fence']
html_theme = 'furo'
html_title = 'Wren 文档'
html_logo = '../../src/wren/assets/wren.svg'
html_favicon = html_logo
html_static_path = ['_static']
html_css_files = ['wren.css']
language = 'zh_CN'
html_search_options = {'dict': str(Path(__file__).with_name('search-terms.txt'))}
exclude_patterns = ['_build']
html_theme_options = {'light_css_variables': {'color-brand-primary': '#116574', 'color-brand-content': '#116574'}}
