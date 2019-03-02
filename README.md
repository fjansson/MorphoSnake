# MorphoSnake

MorphoSnake is a leaf and thallus morphometry program written in Python.

## Installation

### Ubuntu and Debian-like Linux distributions

    sudo apt-get install git python3-numpy python3-scipy python3-matplotlib python3-pip python3-networkx \
    libjpeg-dev libpng-dev libpng12-dev libtiff4-dev libwebp-dev xcftools 

    sudo pip3 install mahotas imread scikit-image    

### Arch Linux

    sudo pacman -S python-pip python-numpy python-scipy python-matplotlib python-networkx
    sudo pip install mahotas imread scikit-image

Morphosnake should work with both networkx 1.x and 2.x.


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

## Contributors

* Fredrik Jansson (fjansson at abo dot fi) - initial implementation
* Catherine Reeb - methodology, measurement definitions, testing
* Pirom Konglerd - measurement methods
* Roque Giordano - updates for NetworkX 2.x



