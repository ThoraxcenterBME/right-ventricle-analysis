#include <framework/mesh.h>

Mesh loadMeshRV(std::istream& in);
void loadRingFromFile(std::string& name, std::vector<Vertex>& vertices);
void mark_excluded(std::string& in, std::vector<Vertex>& vertices); 
void mark_regions(std::string& in, std::vector<Vertex>& vertices); 
char* getOpenFileName(char* szfilter, char* DataPath);
int GetNrOfFiles(std::string fileName);