## 2. Ring statistics analysis

There are several different definitions of a ring in atom system. It is important to make it clear before any further comparison of the data. This is because it determines how the analysis is performed and affects the results significantly. There are some definitions that have been used:

**King’s criteria** - given two adjacent nodes (atoms), a ring is the shortest path from one node (atom) to the other.

**Shortest path criteria** - given a center node (atom) and two of its neighbors, a ring is the shortest path from one neighbor to the other without passing the center node (atom).

**Primitive ring criteria**^[https://doi.org/10.1016/0022-3093(90)90686-G]^[https://doi.org/10.1016/S0927-0256(01)00256-7] - a primitive ring is that the two paths between any node (atom) on the ring and its corresponding prime mid-node (the furthest node on the ring) are both shortest paths. In practice, larger rings can be determined.

**Strong ring criteria**^[https://doi.org/10.1016/0022-3093(91)90145-V] - a ring that can not be decomposed as a sum of any number of smaller rings is a strong ring. Definition of strong ring is  tricky, e.g., it is intuitional that Ring 1 in is a strong ring, and Ring 3 is not, but by flipping point A and B, Ring 2 and Ring 3 can be swapped, even Ring 1 can be regarded as topologically encloses Ring 2 and Ring 3. Hence, it is hard to understand how to identify a strong ring.

In this tool, the Primitive Ring criteria is applied. The analysis needs a parameter to limit the maximum size of the ring, and also limit the memory cost. It is performed in following steps.

### 2.1 Simple ring statistics analysis method

This method uses a simple and straightforward method to identify the primitive rings. The shortest paths list is created, the rings are formed and then checked if any shortcut can be found. If not, they are considered primitive.

There are two problems in earlier developed code:

1. The code gets slower when the ring list gets larger and new primitive rings are kept pushing into the ring list. This is because the ring list has a dynamic size, and allocating memory can be slower when the size is large. Of course, dumping the rings found to a external file can solve this problem, but in this way it is not impossible to avoid repeat analysis of ring.\
To improve the performance, the ring list size is set to a constant and will double its size when full.

2. Many repeat analysis are performed. Usually the repeat analysis of a ring is regarded as a double check of if the ring is really primitive. This is more like a compromise.\
To solve this problem, for each ring found in the code, the atoms on the ring are stored in two ways: unsorted (to keep the topological information), and sorted (to identity each unique ring). When a new ring is found and not checked, the search is first performed to make sure that no repeat analysis is performed.

In the following are the detailed steps of RSA.

#### 2.1.1 Create shortest path array

First, given a center atom, and the already constructed neighbor list, a breadth-first algorithm is used to create a shortest path array about the center atom. It stores all the possible shortest paths begin from the center to the allowed furthest atoms.

However, there are still a few of paths that are not shortest paths. They are the branches of odd rings, which means the end of the path is linked to the second last element of another path. These paths are skipped, waiting for processing in next step.

#### 2.1.2 Construct rings

In this step two types of rings are checked, i.e., even ring and odd ring. For an even ring, if an atom, namely *A*, appears twice or more in the same column (*n*) of the path list array, it is an even ring. For an odd ring, if atom *A* appears in column *n* and *n-1*, and correspondingly an atom *B* appears in the same row but inverse columns, it is an odd ring.

Visibility array contains whether a path is "visible" to another path, i.e., whether a ring can be formed with these two paths. Give *N* as the total number of paths in the Path list array, the visibility array would be a logical type array with size of *N*x*N*. For two paths namely *a* and *b*, the logical element at position (*a*,*b*) determines whether these two paths can form a ring. It has following benefits: 1. From last step we know the odd ring branches should not have further elements, so we just set their visibility false to all other paths. 2. For even rings, all the branches inherited from one of the ring's two branches should not form new rings, thus the corresponding section of the visibility array should be false. 3. It guarantees that every ring is formed by two distinct branches from the center atom.

The visibility array is the part that costs the largest memory, though it is just an array of logical type. When the atom coordinates are large in general, the path list will usually has a large size. Assuming the size is *K*, then the visibility array will have a size of *KxK*. If the memory cost of path list array is on the level of tens of MB, the visibility array can cost hundreds of GB! So this part is the bottle neck of memory at this moment. According to the primitive ring criteria, it is necessary that between a pair of nodes there must be two shortest paths. Therefore, from the shortest path array, the atoms that appear for more than once are considered end nodes of rings. Each ring found is further checked before recognized as a primitive ring.

#### 2.1.3 Remove non-primitive rings

For each ring found in last step, we first check that it is not found already. Because of the sorted index of atoms on the ring, a binary search method can be performed to make search fast. If the max ring size is `M`, and the number of rings is `N`, the complexity of this part is roughly `O(MlogN)`.

If the ring is not analyzed yet, a shortest path search is performed on each node and its mid-node, as defined by Yuan. This is to check if there is shortcut other than the two ring branches between them. If any shortcut is found, the ring will not be pushed to the is removed from the ring list, so that the rings left in the list are all primitive rings at the end of analysis.

The primitive ring search will be performed over all the atoms to find all the primitive rings.

### 2.2 Another method

There is another routine to perform the RSA, the steps are shown below. This method is still in development.

#### 2.2.1 Create shortest path array

The shortest path array (SPA) contains all the possible shortest paths starts from the
center node. Basically, it is an array that contains all the breadth-first search
results. It guarantees that all the rings identified are formed with two
shortest paths.

#### 2.2.2 Modify visibility array

Visibility array (VA) determines whether two shortest paths are ‘visible’ to each other.
At the beginning, most shortest paths are visible to each other, except the path and itself. Thus the VA is a matrix of all `.true.` but the diagonal elements are `.false.`. Given a ring, the VA can be modified in such a way that, the two branches of the ring and all the subsequent branches are 'invisible' to each other, which means all the corresponding elements are set to `.false.`.

VA is useful in two ways:

1. It can avoid repeat identifying of rings around a center atom. Before the analysis about the center atom, this is achieved by modifying the VA. The rings that contains the center atom are found, and the VA is modified according to the rings already found.

2. It can avoid further identifying of rings when a ring is already identified. During the analysis, the branches and subsequent branches are set 'invisible' to each other so that further check is avoided.

In practice, the way the VA is modified is also affected by what kind of rings we are looking for.

#### 2.2.3 Construct rings

If two shortest paths meet at the same node, a ring is identified. Then we
can set the two paths (or their upstream paths) invisible to each other to limit the
identifying of new rings. Ideally, at some point the visibility array will be fulfilled
with `.false.`, which means none of the shortest path is visible to any of them, and
we can stop searching for this center node.

The SPA can be large (~10^4 or ~10^5 shortest paths), I tried to dynamically modify
the SPA and VA so that the size of the arrays can be limited. However, in such a
way step 2 and 3 are merged and it is impossible to initialize the VA with known
rings, i.e., can not avoid repetition.

#### 2.2.4 Remove non-primitive rings

### 2.3 Yuan & Comack's algorithm

#### 2.3.1 Compute reference distance map

First pick several atoms and perform Dijkstra's algorithm to get the distance from every atom to these atoms.

#### 2.3.2 Find rings for each atom

For each atom, first we perform Dijkstra's algorithm to give atoms within certain cutoff a distance number. Then, those atoms are checked to find the *mid-nodes*. For an odd ring, the mid-node should have one neighbor with smaller distance than the atom. For an even ring, it should have at least two neighbors that have smaller distance than it.

When the mid-nodes are found, we follows the decrease of the distance to get the branches from the center node to current node. A ring can be formed for later check.

#### 2.3.3 Check each ring for shortcut

Every pair of atom and its corresponding mid-node on the found ring is checked to make sure there is no short-cut between them. Otherwise they are removed from the list.
