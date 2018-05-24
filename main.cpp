#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <iostream>
#include <igl/triangle_triangle_adjacency.h>
#include <queue>
#include <set>
#include <igl/find.h>
#include <igl/remove_unreferenced.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::MatrixXd N_vertices;
Eigen::MatrixXd N_faces;
Eigen::MatrixXd N_corners;


// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	switch (key)
	{
	case '1':
		viewer.data().set_normals(N_faces);
		return true;
	case '2':
		viewer.data().set_normals(N_vertices);
		return true;
	case '3':
		viewer.data().set_normals(N_corners);
		return true;
	default: break;
	}
	return false;
}

int findFirstZero(Eigen::VectorXi v) {
	for (int i = 0; i < v.rows(); i++)
		if (v[i] == 0)
			return i;
	return -1;
}

void splitMesh(Eigen::MatrixXd V, Eigen::MatrixXi F, int bound = 65000) {

	// Plot the mesh
	std::vector<igl::opengl::glfw::Viewer> viewers;

	// triangle triangle adjacency
	Eigen::MatrixXi TT;
	igl::triangle_triangle_adjacency(F, TT);

	// face discover to record if the face is already into the queue
	//int *FD = new int[F.cols()]{ 0 };
	Eigen::VectorXi FD = Eigen::VectorXi::Zero(F.rows());

	// face queue
	std::queue<int> FQ;

	// if there is still face not being discovered.
	while (!FD.isOnes()) {
		// find the first non-one
		int curIdx = findFirstZero(FD);
		if (curIdx < 0) break;

		// visit it
		FQ.push(curIdx);

		// set of face and vertex for current submesh
		std::set<int> FS;
		std::set<int> VS;

		while (!FQ.empty()) {
			int curFace = FQ.front();
			FD(curFace) = 1;
			FQ.pop();
			//std::cout << "\ndealing " << curFace << " pushing neighbours ";

			FS.insert(curFace);
			for(int i = 0; i < F.cols(); i++)
				VS.insert(F(curFace,i));

			// get three adjacencies and push into the queue aka visit
			for (int i = 0; i < F.cols(); i++) {
				int neighbour = TT(curFace, i);

				// check if already discoved
				if ((neighbour != -1) &&( FD(neighbour) == 0)) {
					FQ.push(neighbour);
					//std::cout << neighbour << "\t";
				}
				//std::cout << "\n";
			}

			// check the size of vertices and faces
			if (FS.size() > bound || VS.size() > bound)
				break;
		}

		// deal with current submesh, turn set of faces to MatrixXd
		Eigen::MatrixXd subV;
		Eigen::MatrixXi subF(FS.size(), F.cols());
		int i = 0;
		for (std::set<int>::iterator it = FS.begin(); it != FS.end(); ++it, i++) {
			subF.row(i) = F.row(*it);
		}
		Eigen::VectorXi UJ;
		igl::remove_unreferenced(V, subF, subV, subF, UJ);


		// write down
		igl::writeOFF("split" + std::to_string(viewers.size()) + ".off", subV, subF);
		std::cout << "split once\n";

		// Plot the mesh
		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(subV, subF);
		viewers.push_back(viewer);

		
	}

	for each (igl::opengl::glfw::Viewer v in viewers)
	{
		v.launch();
	}
}

int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	igl::readOFF( "D:\\Projects\\libigl\\tutorial\\shared\\planexy.off", V, F);

	splitMesh(V, F, 20);

	// Compute per-face normals
	igl::per_face_normals(V, F, N_faces);

	// Compute per-vertex normals
	igl::per_vertex_normals(V, F, N_vertices);

	// Compute per-corner normals, |dihedral angle| > 20 degrees --> crease
	igl::per_corner_normals(V, F, 20, N_corners);

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.callback_key_down = &key_down;
	viewer.data().show_lines = false;
	viewer.data().set_mesh(V, F);
	viewer.data().set_normals(N_faces);
	std::cout <<
		"Press '1' for per-face normals." << std::endl <<
		"Press '2' for per-vertex normals." << std::endl <<
		"Press '3' for per-corner normals." << std::endl;
	viewer.launch();
}