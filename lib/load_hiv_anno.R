# Load HIV Annotations
# 
# Author: gimholte
###############################################################################

HIV_env <- loadFeatures()
lsCategory(HIV_env)
loops <- getFeature(HIV_env, c("loop", "membrane"))
loop_track = ATrack(start = start(loops), end = end(loops), id = getName(loops),
        name = "regions", shape = "box", showFeatureId = TRUE,
        rotation = 45, fontsize = 9)
axis_track = ProteinAxisTrack(addNC = TRUE, littleTicks = TRUE)

gp = getFeature(HIV_env, "protein")
gp_track = ATrack(start = start(gp), end = c(508, 857), id = getName(gp),
        name = "protein", shape = "box", showFeatureId = TRUE,
        fontsize = 14, fontcolor.item = "black", fill = "lightblue")
displayPars(loop_track) = list(fontcolor.item = "black",
        fill = "lightblue")


