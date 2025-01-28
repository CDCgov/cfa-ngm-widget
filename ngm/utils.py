def math_block_format(
    source: str,
    language: str,
    css_class: str,
    options,
    md,
    classes=None,
    id_value="",
    attrs=None,
    **kwargs,
) -> str:
    """Format a math superfence block

    See mkdocs.yml for where this is called. This function was needed
    because ```math blocks are treated as code blocks unless the "math"
    part is specifically recognized and parsed differently by mkdocs.

    Args:
        source (str): content of the block
        language (str): e.g., "math" if the block is ```math
        css_class (str): class of the resulting span
        options (_type_): _description_
        md (_type_): _description_
        classes (_type_, optional): _description_. Defaults to None.
        id_value (str, optional): _description_. Defaults to "".
        attrs (_type_, optional): _description_. Defaults to None.

    Returns:
        str: html to be used in place of markdown
    """
    assert css_class == "arithmatex"
    assert language == "math"
    return f"<span class='{css_class}'> $$ {source} $$ </span>"
