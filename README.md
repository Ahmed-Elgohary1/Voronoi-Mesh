# Voronoi-Mesh Generator

This C++ program generates a Voronoi mesh using Voronoi cells, Delaunay triangulation, and the [Maximal Poisson-Disk Sampling problem](https://github.com/Ahmed-Elgohary1/Mesh-generation). The program takes a set of points as input and generates a Voronoi diagram, which is a partitioning of a plane into regions based on the distance to points in a specific subset of the plane. The Voronoi diagram is then used to create a mesh, which is a collection of vertices, edges, and faces that define a 2D object.

It can make Voronoi mesh for any shape and separate it into domains while decreasing the time complexity and magnitude of memory used
compared to commonly used techniques. 🦄


![image](https://user-images.githubusercontent.com/67281513/163812327-268938fe-f250-46a2-a6cc-026a57ef0fe2.png)


## Installation

To install the program, follow these steps:

1. Clone the repository:
``` git clone https://github.com/Ahmed-Elgohary1/Voronoi-Mesh.git```
2. Build the program using CMake: cmake . && make

