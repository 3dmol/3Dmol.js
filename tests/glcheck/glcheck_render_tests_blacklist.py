# Numbered Tests
# Does not render at all: 
no_render = [
    "test2", 
    "test3", 
    "test5", 
    "test6"
]

# Missing numbers:
missing_numbers = [
    "test35", 
    "test67"
]

#Semi-opaque blob not rendering correctly: 
semi_opaque_surface_not_rendering = [
    "test19", 
    "test24", 
    "test25", 
    "test26", 
    "test27", 
    "test45", 
    "test57", 
    "test58", 
    "test59"
]

#Parts of images not folded/oriented correctly: 
incorectly_folded = [
    "test36", 
    "test38", 
    "test39", 
    "test53"
]

# Text bigger & Colors Deeper: 
text_or_colors_wrong  = [
    "test9", 
    "test69"
]

# Massive size difference: 
size_issues = [
    "test55"
]

# Dark/thicker part of line doesn't render correctly:
line_rendering_issues = [
    "test77"
]





blacklist = [
    *no_render,
    *missing_numbers,
    *semi_opaque_surface_not_rendering,
    *incorectly_folded,
    *text_or_colors_wrong,
    *size_issues,
    *line_rendering_issues
]