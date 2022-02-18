# TODO: determine if these are test fails
questinable_tests = [
    "glmodel_js_setCoordinates"
]

# animations:
animation_failure = [
    "animatelabels",
    "animatereslabels",
    "animateshapes",
    "js_glviewer_js_addModelsAsFrames",
    "js_glviewer_js_zoom",
    "js_glviewer_js_zoomTo",
    "js_glmodel_js_setCoordinates",
    "test90",
    "testvibratearrows",
    "testvibrateboth",
]

# Does not render at all: 
no_render = [
    "embed_zoomto",
    "residiv",
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
surface_not_rendering = [
    "glviewer_addModel",
    "clicksphere",
    "specs_SurfaceStyleSpec",
    "test19", 
    "test24", 
    "test25", 
    "test26", 
    "test27", 
    "test45", 
    "test57", 
    "test58", 
    "test59",
    "test87",
    "test97",
    "test99",
    "testAddModelPQR",
    "testallselsurf",
    "testSurfaceVolData",
    "testsymmsurf",
    "test_volumetric_seldist",
    "testCIFsurf"
]

geometry_not_rendering = [
    "refreshvibrate", #little arrow missing
    "testbackalpha",
    "testmulticif",
]

extra_geomety = [
    "testVolumetricRenderCorners",
    "testVolumetricRenderRod",
]

#Parts of images not folded/oriented correctly: 
incorectly_folded = [
    "test36", 
    "test38",
    "test39", 
    "test53"
]

# Massive size difference: 
size_issues = [
    "js_glviewer_js_addLabel",
    "test55",
    #"testmulticif"
]

# Dark/thicker part of line doesn't render correctly:
line_rendering_issues = [
    "test77",
    "test93",
    "test91",
    "test94"
]

blacklist = [
    *no_render,
    *missing_numbers,
    *surface_not_rendering,
    *incorectly_folded,
    *extra_geomety,
    *size_issues,
    *line_rendering_issues,
    *animation_failure,
]


print(f"the blacklist currently contains {len(blacklist)} tests")
