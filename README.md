# MorphoSnake

MorphoSnake is a leaf and thallus morphometry program written in Python.

Authors: Fredrik Jansson (fjansson at abo dot fi),
	 Pirom Konglerd,
         Catherine Reeb



## Installation

### Ubuntu and Debian-like Linux distributions

    sudo apt-get install git python3-numpy python3-scipy python3-matplotlib python3-pip python3-networkx \
    libjpeg-dev libpng-dev libpng12-dev libtiff4-dev libwebp-dev xcftools 

    sudo pip3 install mahotas imread scikit-image    

### Arch Linux

    sudo pacman -S python-pip python-numpy python-scipy python-matplotlib python-networkx
    sudo pip install mahotas imread scikit-image

Please note: MorphoSnake currently depends on networkx version 1.x, while some distributions already ship 2.0.
Until MorphoSnake is updated, networkx 1.x can be installed locally, e.g. in a virtual environment, with pip3 install 'networkx<2.0'.


## Usage and Measurements

MorphoSnake operates on images of leaves or other branched structures, and is used to extract features from the structure.

The image is skeletonized using the Zhang-Suen skeletonization method, as implemented in [scikit image](http://scikit-image.org/docs/dev/api/skimage.morphology.html#skimage.morphology.skeletonize)

From the skeleton, a graph is constructed. The graph has vertices at the branching points of the skeleton.
Since the image skeleonization method is not perfect, some spurious branches may appear. Also, if different branches overlap,
the generated skeleton can contain loops. The loops and the spurious branches can be removed manually, by deleting
nodes or vertices from the graph.

Each edge and node is in the graph assigned an index according to the Strahler hierarchy: terminal branches ("tips") are assigned number 1. The parent branch of two branches of level n is assigned level n+1. The parent of branches of different level is assigned the same level as the child branch of the highest level.

For the edges and nodes in the graph, the following measurements are collected.

#### Node
* `level`    - distance in the graph to the nearest leaf. A leaf has level 1.
* `Strahler` - Strahler order. 
* `dia`      - diameter

#### Edge
* `level`          - distance in the graph to the nearest leaf. A leaf has level 1.
* `Strahler`       - Strahler order. 
* `branchlength`   - length of the branch along the skeleton
* `branchlength_e` - Euclidian (straight-line) distance between the end points  
* `alpha`          - angle from parent branch, measured at the node's disk
* `alpha_e`        - angle from parent branch, measured using the straight graph edges
* `W_mean`         - mean width of the branch
* `W_max`          - max width of the branch
* `W_min`          - min width of the branch
* `W_std`          - standard deviation of the width of the branch



## Image resolution data base

In order to report measurement results in physical units, the program needs to know the resolution of the input image.
The resolution, in pixels per mm, and optionally other properties, are stored in a JSON file, `leaves.json`, in the same directory as the MorphoSnake script.

    "imagename": {
        "px_mm": 430.0
    },



When the program is run on a leaf image, it stores the location of the root node in the same data base. This lets the
program remember the selected root for each image beteen runs of the program.


### Keyboard commands

* 'd' - delete node closest to pointer
* 'x' - delete edge closest to pointer
* 'u' - undo last action
* 'r' - rebuild the graph - restoring all manually deleted nodes and edges
* click a node to select it as the root

### License

MorphoSnake is released under the GNU General Public License v3.0,

see the file LICENSE for details.

## Updated version

Installation: (Assuming you have all the requirements for the old morphosnake version) 

to run the updated version of morphosnake you have to update networkx to version 2.0 via:

sudo pip3 install networkx --upgrade

The library OpenCV (cv2) is also needed, you can install it using: 

sudo pip3 install opencv-python

To charge an image, make sure no graph file (.pkl)created by an older version is available. please delete it and try again
Keyboard commands:

Click on a node to select it, then use one of the following commands:

	't' to change the root to the selected node 
	'a' to mark the selected node as fertile
	'd' to delete the selected node
	'b' to modify the diameter of the node (just place the mouse where the new diameter is wanted and press the button)
	
Without any click, press on the keyboard:

	'n' to add a node to the skeleton pixel closest to the pointer
	'x' to delete edge closest to the pointer
	'u' to undo last modification
	'r' to rebuild the graph - restoring all manually deleted nodes and edges
	'alt+e' to hide the changes graphically
	




