# tpf2-make-walschaerts-ani
A script tool to generate .ani file for rod animation.

# How to use
1. Perpare models.
2. Set parameters in config.lua.
3. Execute MakeWalschaertsAni.exe
4. Set "animations" block to your .mdl like this:

```
    {
        animations = {
            drive = {
                params = {
                    id = "vehicle/train/jgr_8620/connecting_rod.ani",
                },
                type = "FILE_REF",
            },
        },
        materials = { "vehicle/train/jgr_8620/base.mtl", },
        mesh = "vehicle/train/jgr_8620/loco_connecting_rod_lod0.msh",
        name = "loco_connecting_rod_right",
        transf = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, -1, 0.49499997496605, 1, },
    },
    {
        animations = {
            drive = {
                params = {
                    id = "vehicle/train/jgr_8620/connecting_rod.ani",
                },
                type = "FILE_REF",
            },
        },
        materials = { "vehicle/train/jgr_8620/base.mtl", },
        mesh = "vehicle/train/jgr_8620/loco_connecting_rod_lod0.msh",
        name = "loco_connecting_rod_left",
        transf = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0.49499997496605, 1, },
    },
```
