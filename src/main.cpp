/*
 * @file main.cpp
 * @brief This code loads the .OBJ mesh, passes vertices, normals,
 * eye normal and light position to shader
 *
 * @author Vishwesh Vhavle
 * @date September 20, 2023
 */

#include "utils.h"
#include <vector>

using namespace std;

int screen_width = 800, screen_height = 450;
vector<int> voxel_grid;
vector<int> distance_field;

GLint vModel_uniform, vView_uniform, vProjection_uniform;
GLint lpos_world_uniform, eye_normal_uniform;
glm::mat4 modelT, viewT, projectionT;
GLint vColor_uniform;

double oldX, oldY, currentX, currentY;
bool isDragging = false;

void rasterizer();
void createMeshObject(unsigned int &, unsigned int &);
void createQuadMeshObjectWithoutNormals(unsigned int &, unsigned int &);
void createVoxelizedMeshObject(unsigned int &, unsigned int &, int);
void createVoxelizedQuadMeshObject(unsigned int &program, unsigned int &shape_VAO, int grid_size);

void setupModelTransformation(unsigned int &);
void setupViewTransformation(unsigned int &);
void setupProjectionTransformation(unsigned int &);
glm::vec3 getTrackBallVector(double x, double y);

GLuint nVertices;
float scale = 0.05; // Change the scale of the model as needed

struct vert{
    glm::vec3 pos;
    glm::vec3 n;
    glm::vec2 uv;
};

void vert_stats(const vector<glm::vec3>& vertices) {
    float min_x = vertices[0].x;
    float min_y = vertices[0].y;
    float min_z = vertices[0].z;
    float max_x = vertices[0].x;
    float max_y = vertices[0].y;
    float max_z = vertices[0].z;

    for (const glm::vec3& v : vertices) {
        if (v.x < min_x) {
            min_x = v.x;
        }
        if (v.y < min_y) {
            min_y = v.y;
        }
        if (v.z < min_z) {
            min_z = v.z;
        }
        if (v.x > max_x) {
            max_x = v.x;
        }
        if (v.y > max_y) {
            max_y = v.y;
        }
        if (v.z > max_z) {
            max_z = v.z;
        }
    }

    // cout << "min: " << min_x << " " << min_y << " " << min_z << endl;
    // cout << "max: " << max_x << " " << max_y << " " << max_z << endl;

    exit(0);
}

void normalize_vertices(vector<glm::vec3>& vertices) {
    // send all vertices to be between 0 and 1

    // vector<glm::vec3> vertices = vertices_p;
    float min_x = vertices[0].x;
    float min_y = vertices[0].y;
    float min_z = vertices[0].z;
    float max_x = vertices[0].x;
    float max_y = vertices[0].y;
    float max_z = vertices[0].z;

    // calculate min and max
    for (glm::vec3& v : vertices) {
        if (v.x < min_x) {
            min_x = v.x;
        }
        if (v.y < min_y) {
            min_y = v.y;
        }
        if (v.z < min_z) {
            min_z = v.z;
        }
        if (v.x > max_x) {
            max_x = v.x;
        }
        if (v.y > max_y) {
            max_y = v.y;
        }
        if (v.z > max_z) {
            max_z = v.z;
        }
    }

    // normalize
    for (glm::vec3& v : vertices) {
        // cout << v.x << " " << v.y << " " << v.z << endl;
        v.x = (v.x - min_x) / (max_x - min_x);
        v.y = (v.y - min_y) / (max_y - min_y);
        v.z = (v.z - min_z) / (max_z - min_z);

        // cout << v.x << " " << v.y << " " << v.z << endl << endl;
    }

    // cout << "min: " << min_x << " " << min_y << " " << min_z << endl;
    // cout << "max: " << max_x << " " << max_y << " " << max_z << endl;
}

void initiate_voxel_grid(int gridSize) {
    cout << "initializing voxel grid of size " << gridSize << 'x' << gridSize << 'x' << gridSize << " = " << gridSize * gridSize * gridSize << std::endl;

    voxel_grid.resize(gridSize * gridSize * gridSize);
    for (int i = 0; i < gridSize * gridSize * gridSize; i++) {
        voxel_grid[i] = 0;
    }
}

void voxelization(const vector<glm::vec3>& vertices, int gridSize) {
    // for (const glm::vec3& v : vertices) {
    //     int x = (v.x + 0.5f) * gridSize;
    //     int y = (v.y + 0.5f) * gridSize;
    //     int z = (v.z + 0.5f) * gridSize;
    //     cout << x << " " << y << " " << z << endl;
    //     voxel_grid[x + y * gridSize + z * gridSize * gridSize] = 1;
    // }

    for (const glm::vec3& v : vertices) {
        int x = (v.x) * gridSize;
        int y = (v.y) * gridSize;
        int z = (v.z) * gridSize;

        if (v.z > 0.5) {
            // cout << v.z << " " << x << " " << y << " " << z << endl;
        }
        // cout << gridSize * gridSize * gridSize << " " << x + y * gridSize + z * gridSize * gridSize << endl;
        voxel_grid[x + y * gridSize + z * gridSize * gridSize] = 1;
    }

    for (int k = 0; k < gridSize; k++) {
        for (int j = 0; j < gridSize; j++){
            for (int i = 0; i < gridSize; i++) {
                cout << voxel_grid[k * gridSize * gridSize + j * gridSize + i] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

void create_distance_field(int gridSize) {
    cout << "creating (squared) distance field" << endl;
    distance_field.resize(gridSize * gridSize * gridSize);
    for (int i = 0; i < gridSize * gridSize * gridSize; i++) {
        distance_field[i] = 0;
    }

    for (int k = 0; k < gridSize; k++) {
        
        for (int j = 0; j < gridSize; j++){
            
            for (int i = 0; i < gridSize; i++) {
                
                if (voxel_grid[k * gridSize * gridSize + j * gridSize + i] == 1) {
                    distance_field[k * gridSize * gridSize + j * gridSize + i] = 0;
                }
                else {
                    int min_dist = gridSize * gridSize * gridSize;
                    for (int k2 = 0; k2 < gridSize; k2++) {
                        for (int j2 = 0; j2 < gridSize; j2++){
                            for (int i2 = 0; i2 < gridSize; i2++) {
                                if (voxel_grid[k2 * gridSize * gridSize + j2 * gridSize + i2] == 1) {
                                    int dist = (k - k2) * (k - k2) + (j - j2) * (j - j2) + (i - i2) * (i - i2);
                                    if (dist < min_dist) {
                                        min_dist = dist;
                                    }
                                }
                            }
                        }
                    }
                    distance_field[k * gridSize * gridSize + j * gridSize + i] = min_dist;
                }

            }
            cout << "k = " << k << " j = " << j << endl;
        }
    }

    for (int k = 0; k < gridSize; k++) {
        for (int j = 0; j < gridSize; j++){
            for (int i = 0; i < gridSize; i++) {
                if (distance_field[k * gridSize * gridSize + j * gridSize + i] < 10) {
                    cout << "  " << distance_field[k * gridSize * gridSize + j * gridSize + i] << " ";
                }
                else if (distance_field[k * gridSize * gridSize + j * gridSize + i] < 100) {
                    cout << " " << distance_field[k * gridSize * gridSize + j * gridSize + i] << " ";
                }
                else {
                    cout << distance_field[k * gridSize * gridSize + j * gridSize + i] << " ";
                }
            }
            cout << endl;
        }
        cout << endl;
    }
}

void get_cube_vertex_array(int i, int j, int k, vert* vertices, int gridSize) {
    vertices[ 0].pos = glm::vec3(i,   j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 1].pos = glm::vec3(i+1, j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 2].pos = glm::vec3(i+1, j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 3].pos = glm::vec3(i+1, j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 4].pos = glm::vec3(i,   j  , k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 5].pos = glm::vec3(i,   j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);

    vertices[ 6].pos = glm::vec3(i,   j,   k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 7].pos = glm::vec3(i+1, j,   k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 8].pos = glm::vec3(i+1, j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[ 9].pos = glm::vec3(i+1, j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[10].pos = glm::vec3(i,   j  , k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[11].pos = glm::vec3(i,   j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);

    vertices[12].pos = glm::vec3(i,   j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[13].pos = glm::vec3(i+1, j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[14].pos = glm::vec3(i+1, j,   k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[15].pos = glm::vec3(i+1, j,   k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[16].pos = glm::vec3(i,   j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[17].pos = glm::vec3(i,   j,   k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);

    vertices[18].pos = glm::vec3(i,   j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[19].pos = glm::vec3(i+1, j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[20].pos = glm::vec3(i+1, j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[21].pos = glm::vec3(i+1, j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[22].pos = glm::vec3(i,   j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[23].pos = glm::vec3(i,   j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);

    vertices[24].pos = glm::vec3(i,   j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[25].pos = glm::vec3(i,   j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[26].pos = glm::vec3(i,   j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[27].pos = glm::vec3(i,   j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[28].pos = glm::vec3(i,   j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[29].pos = glm::vec3(i,   j,   k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);

    vertices[30].pos = glm::vec3(i+1, j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[31].pos = glm::vec3(i+1, j+1, k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[32].pos = glm::vec3(i+1, j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[33].pos = glm::vec3(i+1, j+1, k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[34].pos = glm::vec3(i+1, j,   k  ) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);
    vertices[35].pos = glm::vec3(i+1, j,   k+1) - glm::vec3(0.5*gridSize, 0.5*gridSize, 0.5*gridSize);

    for (int i = 0; i < 36; i+=3) {
        glm::vec3 v1 = glm::vec3(vertices[(i + 1)].pos.x - vertices[i].pos.x, vertices[(i + 1)].pos.y - vertices[i].pos.y, vertices[(i + 1)].pos.z - vertices[i].pos.z);
        glm::vec3 v2 = glm::vec3(vertices[(i + 2)].pos.x - vertices[(i + 1)].pos.x, vertices[(i + 2)].pos.y - vertices[(i + 1)].pos.y, vertices[(i + 2)].pos.z - vertices[(i + 1)].pos.z);

        // using the first two vectors of a triangle, we can calculate the normal
        glm::vec3 n = norm(v1, v2);

        if (cross(v1, v2).x < 0) {
            n = -n;
        }
        if (cross(v1, v2).y < 0) {
            n = -n;
        }
        if (cross(v1, v2).z < 0) {
            n = -n;
        }


        vertices[i].n = n;
        vertices[i+1].n = n;
        vertices[i+2].n = n;
    }
}



vector<vert> create_vertices_from_voxel_grid(int gridSize) {
    vector<vert> vertices;
    for (int k = 0; k < gridSize; k++) {
        for (int j = 0; j < gridSize; j++){
            for (int i = 0; i < gridSize; i++) {
                // cout << gridSize * gridSize * gridSize << " " << i * gridSize * gridSize + j * gridSize + k << " ";
                if (voxel_grid[k * gridSize * gridSize + j * gridSize + i] == 1) {
                    // cout << " " << "<--------";
                    vert cube_v[36];
                    get_cube_vertex_array(i, j, k, cube_v, gridSize);
                    for (int idx = 0; idx < 36; idx++) {
                        vertices.push_back(cube_v[idx]);
                    }
                // cout << endl;
                }
            }
        }
    }
    return vertices;
}


int main()
{
    int grid_size = 50;
    initiate_voxel_grid(grid_size);

    // Setup window
    GLFWwindow *window = setupWindow(screen_width, screen_height);
    ImGuiIO &io = ImGui::GetIO(); // Create IO object
    ImVec4 clearColor = ImVec4(0.81f, 0.78f, 0.81f, 1.00f);

    unsigned int shaderProgram = createProgram("./shaders/vshader.vs", "./shaders/fshader.fs");
    // unsigned int shaderProgram = createProgram("./shaders_gouraud/vshader.vs", "./shaders_gouraud/fshader.fs");

    // Get handle to light position variable in shader
    lpos_world_uniform = glGetUniformLocation(shaderProgram, "lpos_world");
    if (lpos_world_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: lpos_world\n");
    }

    // Get handle to eye normal variable in shader
    eye_normal_uniform = glGetUniformLocation(shaderProgram, "eye_normal");
    if (eye_normal_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: eye_normal. Specular Lighting Switched Off.\n");
    }

    // Get handle to color variable in shader
    vColor_uniform = glGetUniformLocation(shaderProgram, "vColor");
    if (vColor_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: vColor\n");
        exit(0);
    }

    glUseProgram(shaderProgram);

    unsigned int VAO;
    glGenVertexArrays(1, &VAO);

    setupModelTransformation(shaderProgram);
    setupViewTransformation(shaderProgram);
    setupProjectionTransformation(shaderProgram);

    // createMeshObject(shaderProgram, VAO);
    // createQuadMeshObjectWithoutNormals(shaderProgram, VAO);
    createVoxelizedMeshObject(shaderProgram, VAO, grid_size);

    oldX = oldY = currentX = currentY = 0.0;
    int prevLeftButtonState = GLFW_RELEASE;
    // int mode = 1; // 0 for manual drag, 1 for auto rotate
    int mode = 0;

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        double current_seconds = glfwGetTime();

        // Get current mouse position
        int leftButtonState = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
        double x, y;
        glfwGetCursorPos(window, &x, &y);
        if (leftButtonState == GLFW_PRESS && prevLeftButtonState == GLFW_RELEASE)
        {
            isDragging = true;
            currentX = oldX = x;
            currentY = oldY = y;
        }
        else if (leftButtonState == GLFW_PRESS && prevLeftButtonState == GLFW_PRESS)
        {
            currentX = x;
            currentY = y;
        }
        else if (leftButtonState == GLFW_RELEASE && prevLeftButtonState == GLFW_PRESS)
        {
            isDragging = false;
        }

        // Rotate based on mouse drag movement
        float angle = 0.0f;
        prevLeftButtonState = leftButtonState;
        if (mode == 0)
        {
            if (isDragging && (currentX != oldX || currentY != oldY))
            {
                // Drag rotation
                glm::vec3 va = getTrackBallVector(oldX, oldY);
                glm::vec3 vb = getTrackBallVector(currentX, currentY);

                angle = acos(std::min(1.0f, glm::dot(va, vb)));
                glm::vec3 axis_in_camera_coord = glm::cross(va, vb);
                glm::mat3 camera2object = glm::inverse(glm::mat3(viewT * modelT));
                glm::vec3 axis_in_object_coord = camera2object * axis_in_camera_coord;
                modelT = glm::rotate(modelT, angle, axis_in_object_coord);
                glUniformMatrix4fv(vModel_uniform, 1, GL_FALSE, glm::value_ptr(modelT));

                oldX = currentX;
                oldY = currentY;
            }
        }
        else if (mode == 1)
        {
            // Autorotation
            angle += 0.05;
            modelT = glm::rotate(modelT, angle, glm::vec3(0.0, 1.0, 0.0));
            glUniformMatrix4fv(vModel_uniform, 1, GL_FALSE, glm::value_ptr(modelT));
        }

        // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        glUseProgram(shaderProgram);

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clearColor.x, clearColor.y, clearColor.z, clearColor.w);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glBindVertexArray(VAO);

        glUniform3f(lpos_world_uniform, -50.0, 500.0, 30.0);
        glUniform3f(eye_normal_uniform, -0.0, 0.0, -40.0);
        glUniform3f(vColor_uniform, 0.1, 0.8, 0.85);

        glDrawArrays(GL_TRIANGLES, 0, nVertices);
        // glDrawArrays(GL_QUADS, 0, nVertices);

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // Cleanup
    cleanup(window);
}

void setupModelTransformation(unsigned int &program)
{
    // Modelling transformations (Model -> World coordinates)
    modelT = glm::translate(glm::mat4(1.0f), glm::vec3(0.0, 0.0, 0.0)); // Model coordinates are the world coordinates

    // Pass on the modelling matrix to the vertex shader
    glUseProgram(program);
    vModel_uniform = glGetUniformLocation(program, "vModel");
    if (vModel_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: vModel\n");
        exit(0);
    }
    glUniformMatrix4fv(vModel_uniform, 1, GL_FALSE, glm::value_ptr(modelT));
}

void setupViewTransformation(unsigned int &program)
{
    // Viewing transformations (World -> Camera coordinates
    // Camera at (40, 20, 40)  in a right handed coordinate system
    viewT = glm::lookAt(glm::vec3(60.0, 60.0, 60.0), glm::vec3(30.0, 30.0, 30.0), glm::vec3(0.0, 1.0, 0.0));

    // Pass-on the viewing matrix to the vertex shader
    glUseProgram(program);
    vView_uniform = glGetUniformLocation(program, "vView");
    if (vView_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: vView\n");
        exit(0);
    }
    glUniformMatrix4fv(vView_uniform, 1, GL_FALSE, glm::value_ptr(viewT));
}

void setupProjectionTransformation(unsigned int &program)
{
    // Projection transformation
    projectionT = glm::perspective(45.0f, (GLfloat)screen_width / (GLfloat)screen_height, 0.1f, 1000.0f);

    // Pass on the projection matrix to the vertex shader
    glUseProgram(program);
    vProjection_uniform = glGetUniformLocation(program, "vProjection");
    if (vProjection_uniform == -1)
    {
        fprintf(stderr, "Could not bind location: vProjection\n");
        exit(0);
    }
    glUniformMatrix4fv(vProjection_uniform, 1, GL_FALSE, glm::value_ptr(projectionT));
}

glm::vec3 getTrackBallVector(double x, double y)
{
    glm::vec3 p = glm::vec3(2.0 * x / screen_width - 1.0, 2.0 * y / screen_height - 1.0, 0.0); // Normalize to [-1, +1]
    p.y = -p.y;                                                                                // Invert Y since screen coordinate and OpenGL coordinates have different Y directions.

    float mag2 = p.x * p.x + p.y * p.y;
    if (mag2 <= 1.0f)
        p.z = sqrtf(1.0f - mag2);
    else
        p = glm::normalize(p); // Nearest point, close to the sides of the trackball
    return p;
}

void createQuadMeshObjectWithoutNormals(unsigned int &program, unsigned int &shape_VAO)
{
    vector<int> vertex_indices, uv_indices;
    vector<glm::vec3> temp_vertices;
    vector<glm::vec2> temp_uvs;
    vector<glm::vec3> temp_normals;
    int ct = 0;

    scale = 4; // Change Scale of the model as needed
    // smpl.obj takes ~ 30 as scale, buddha and bunny take ~ 0.05 as scale
    // FILE * file = fopen("src/smpl.obj", "r");
    FILE *file = fopen("src/Cobblestones3.obj", "r");
    if (file == NULL)
        printf("File not found\n");

    while (true)
    {

        char head[128];
        // read the first word of the line
        int res = fscanf(file, "%s", head);
        if (res == EOF)
            break; // EOF = End Of File. Quit the loop.

        if (strcmp(head, "v") == 0)
        {
            glm::vec3 vertex;
            fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
            // printf("%f %f %f\n", vertex.x, vertex.y, vertex.z );
            temp_vertices.push_back(vertex);
            ct++;
        }

        else if (strcmp(head, "vt") == 0)
        {
            glm::vec2 uv;
            fscanf(file, "%f %f\n", &uv.x, &uv.y);
            // printf("%f %f\n", uv.x, uv.y );
            temp_uvs.push_back(uv);
        }

        else if (strcmp(head, "vn") == 0)
        {
            glm::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
            // printf("%f %f %f\n", normal.x, normal.y, normal.z );
            temp_normals.push_back(normal);
        }

        else if (strcmp(head, "f") == 0)
        {
            string vertex1, vertex2, vertex3;
            int vertexIndex[4], uvIndex[4];
            int matches = fscanf(file, "%d/%d %d/%d %d/%d %d/%d\n", &vertexIndex[0], &uvIndex[0], &vertexIndex[1], &uvIndex[1], &vertexIndex[2], &uvIndex[2], &vertexIndex[3], &uvIndex[3]);
            if (matches != 8)
            {
                printf("OBJ File may not contain texture coordinates\n");
            }
            vertex_indices.push_back(vertexIndex[0]);
            vertex_indices.push_back(vertexIndex[1]);
            vertex_indices.push_back(vertexIndex[2]);
            vertex_indices.push_back(vertexIndex[0]);
            vertex_indices.push_back(vertexIndex[2]);
            vertex_indices.push_back(vertexIndex[3]);
            uv_indices.push_back(uvIndex[0]);
            uv_indices.push_back(uvIndex[1]);
            uv_indices.push_back(uvIndex[2]);
            uv_indices.push_back(uvIndex[0]);
            uv_indices.push_back(uvIndex[2]);
            uv_indices.push_back(uvIndex[3]);
        }
    }
    fclose(file);

    glUseProgram(program);

    // Bind shader variables
    int vVertex_attrib = glGetAttribLocation(program, "vVertex");
    if (vVertex_attrib == -1)
    {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    int vNormal_attrib = glGetAttribLocation(program, "vNormal");
    if (vNormal_attrib == -1)
    {
        std::cout << "Could not bind location: vNormal\n";
    }

    GLfloat *shape_vertices = new GLfloat[vertex_indices.size() * 3];
    GLfloat *vertex_normals = new GLfloat[vertex_indices.size() * 3];

    nVertices = vertex_indices.size() * 3;

    // Ordered vertex array to be passed to the shader to draw the mesh
    for (int i = 0; i < vertex_indices.size(); i++)
    {
        int vertexIndex = vertex_indices[i];
        shape_vertices[i * 3] = temp_vertices[vertexIndex - 1][0] * scale; // x-coordinate
        shape_vertices[i * 3 + 1] = temp_vertices[vertexIndex - 1][1] * scale; // y-coordinate
        shape_vertices[i * 3 + 2] = temp_vertices[vertexIndex - 1][2] * scale; // z-coordinate
    }

    // generate normals for the triangle mesh
    for (int i = 0; i < vertex_indices.size(); i += 3)
    {
        glm::vec3 v1 = glm::vec3(shape_vertices[(i + 1) * 3] - shape_vertices[i * 3], shape_vertices[(i + 1) * 3 + 1] - shape_vertices[i * 3 + 1], shape_vertices[(i + 1) * 3 + 2] - shape_vertices[i * 3 + 2]);
        glm::vec3 v2 = glm::vec3(shape_vertices[(i + 2) * 3] - shape_vertices[(i + 1) * 3], shape_vertices[(i + 2) * 3 + 1] - shape_vertices[(i + 1) * 3 + 1], shape_vertices[(i + 2) * 3 + 2] - shape_vertices[(i + 1) * 3 + 2]);

        // using the first two vectors of a triangle, we can calculate the normal
        glm::vec3 n = cross(v1, v2);
        n = normalize(n);


        vertex_normals[i * 3] = n.x;
        vertex_normals[(i + 1) * 3] = n.x;
        vertex_normals[(i + 2) * 3] = n.x;
        vertex_normals[i * 3 + 1] = n.y;
        vertex_normals[(i + 1) * 3 + 1] = n.y;
        vertex_normals[(i + 2) * 3 + 1] = n.y;
        vertex_normals[i * 3 + 2] = n.z;
        vertex_normals[(i + 1) * 3 + 2] = n.z;
        vertex_normals[(i + 2) * 3 + 2] = n.z;
    }

    for (int i = 0; i < vertex_indices.size(); i += 3)
    {

        glm::vec3 a = glm::vec3(shape_vertices[i * 3], shape_vertices[i * 3 + 1], shape_vertices[i * 3 + 2]);
        glm::vec3 b = glm::vec3(shape_vertices[(i + 1) * 3], shape_vertices[(i + 1) * 3 + 1], shape_vertices[(i + 1) * 3 + 2]);
        glm::vec3 c = glm::vec3(shape_vertices[(i + 2) * 3], shape_vertices[(i + 2) * 3 + 1], shape_vertices[(i + 2) * 3 + 2]);
        glm::vec3 n = glm::vec3(vertex_normals[i * 3], vertex_normals[i * 3 + 1], vertex_normals[i * 3 + 2]);

        // std::cout << a.x << " " << a.y << " " << a.z << std::endl;
        // std::cout << b.x << " " << b.y << " " << b.z << std::endl;
        // std::cout << c.x << " " << c.y << " " << c.z << std::endl;
        // std::cout << n.x << " " << n.y << " " << n.z << std::endl << std::endl;
    }

    // Generate VAO object
    glGenVertexArrays(1, &shape_VAO);
    glBindVertexArray(shape_VAO);

    // Create VBOs for the VAO
    GLuint vertex_VBO; // Vertex Buffer
    glGenBuffers(1, &vertex_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertex_indices.size() * 3, shape_vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vVertex_attrib);
    glVertexAttribPointer(vVertex_attrib, 3, GL_FLOAT, GL_FALSE, 00, 0);
    delete[] shape_vertices;

    GLuint normal_VBO; // Normal Buffer
    glGenBuffers(1, &normal_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, normal_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertex_indices.size() * 3, vertex_normals, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vNormal_attrib);
    glVertexAttribPointer(vNormal_attrib, 3, GL_FLOAT, GL_FALSE, 0, 0);
    delete[] vertex_normals;

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}




void createMeshObject(unsigned int &program, unsigned int &shape_VAO)
{
    vector<int> vertex_indices, uv_indices, normal_indices;
    vector<glm::vec3> temp_vertices;
    vector<glm::vec2> temp_uvs;
    vector<glm::vec3> temp_normals;
    int ct = 0;

    scale = 30; // Change Scale of the model as needed
    // scale = 0.05;
    // smpl.obj takes ~ 30 as scale, buddha and bunny take ~ 0.05 as scale
    FILE *file = fopen("src/smpl.obj", "r");
    // FILE *file = fopen("src/bunny.obj", "r");
    // FILE * file = fopen("../Cobblestones3/Files/untitled.obj", "r");
    if (file == NULL)
        printf("File not found\n");

    while (true)
    {

        char head[128];
        // read the first word of the line
        int res = fscanf(file, "%s", head);
        if (res == EOF)
            break; // EOF = End Of File. Quit the loop.

        if (strcmp(head, "v") == 0)
        {
            glm::vec3 vertex;
            fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
            // printf("%f %f %f\n", vertex.x, vertex.y, vertex.z );
            temp_vertices.push_back(vertex);
            ct++;
        }

        else if (strcmp(head, "vt") == 0)
        {
            glm::vec2 uv;
            fscanf(file, "%f %f\n", &uv.x, &uv.y);
            // printf("%f %f\n", uv.x, uv.y );
            temp_uvs.push_back(uv);
        }

        else if (strcmp(head, "vn") == 0)
        {
            glm::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
            // printf("%f %f %f\n", normal.x, normal.y, normal.z );
            temp_normals.push_back(normal);
        }

        else if (strcmp(head, "f") == 0)
        {
            string vertex1, vertex2, vertex3;
            int vertexIndex[3], uvIndex[3], normalIndex[3];
            int matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2]);
            if (matches != 9)
            {
                printf("OBJ File may not contain texture coordinates or normal coordinates\n");
            }
            vertex_indices.push_back(vertexIndex[0]);
            vertex_indices.push_back(vertexIndex[1]);
            vertex_indices.push_back(vertexIndex[2]);
            uv_indices.push_back(uvIndex[0]);
            uv_indices.push_back(uvIndex[1]);
            uv_indices.push_back(uvIndex[2]);
            normal_indices.push_back(normalIndex[0]);
            normal_indices.push_back(normalIndex[1]);
            normal_indices.push_back(normalIndex[2]);
        }
    }
    fclose(file);

    glUseProgram(program);

    // Bind shader variables
    int vVertex_attrib = glGetAttribLocation(program, "vVertex");
    if (vVertex_attrib == -1)
    {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    int vNormal_attrib = glGetAttribLocation(program, "vNormal");
    if (vNormal_attrib == -1)
    {
        std::cout << "Could not bind location: vNormal\n";
    }

    GLfloat *shape_vertices = new GLfloat[vertex_indices.size() * 3];
    GLfloat *vertex_normals = new GLfloat[normal_indices.size() * 3];

    nVertices = vertex_indices.size() * 3;

    for (int i = 0; i < vertex_indices.size(); i++)
    {
        int vertexIndex = vertex_indices[i];
        shape_vertices[i * 3] = temp_vertices[vertexIndex - 1][0] * scale;
        shape_vertices[i * 3 + 1] = temp_vertices[vertexIndex - 1][1] * scale;
        shape_vertices[i * 3 + 2] = temp_vertices[vertexIndex - 1][2] * scale;
    }

    // generated normals for the triangle mesh
    for (int i = 0; i < vertex_indices.size(); i += 3)
    {
        glm::vec3 v1 = glm::vec3(shape_vertices[(i + 1) * 3] - shape_vertices[i * 3], shape_vertices[(i + 1) * 3 + 1] - shape_vertices[i * 3 + 1], shape_vertices[(i + 1) * 3 + 2] - shape_vertices[i * 3 + 2]);
        glm::vec3 v2 = glm::vec3(shape_vertices[(i + 2) * 3] - shape_vertices[(i + 1) * 3], shape_vertices[(i + 2) * 3 + 1] - shape_vertices[(i + 1) * 3 + 1], shape_vertices[(i + 2) * 3 + 2] - shape_vertices[(i + 1) * 3 + 2]);

        glm::vec3 n = norm(v1, v2);

        vertex_normals[i * 3] = n.x;
        vertex_normals[(i + 1) * 3] = n.x;
        vertex_normals[(i + 2) * 3] = n.x;
        vertex_normals[i * 3 + 1] = n.y;
        vertex_normals[(i + 1) * 3 + 1] = n.y;
        vertex_normals[(i + 2) * 3 + 1] = n.y;
        vertex_normals[i * 3 + 2] = n.z;
        vertex_normals[(i + 1) * 3 + 2] = n.z;
        vertex_normals[(i + 2) * 3 + 2] = n.z;
    }

    for (int i = 0; i < vertex_indices.size(); i += 3)
    {

        glm::vec3 a = glm::vec3(shape_vertices[i * 3], shape_vertices[i * 3 + 1], shape_vertices[i * 3 + 2]);
        glm::vec3 b = glm::vec3(shape_vertices[(i + 1) * 3], shape_vertices[(i + 1) * 3 + 1], shape_vertices[(i + 1) * 3 + 2]);
        glm::vec3 c = glm::vec3(shape_vertices[(i + 2) * 3], shape_vertices[(i + 2) * 3 + 1], shape_vertices[(i + 2) * 3 + 2]);
        glm::vec3 n = glm::vec3(vertex_normals[i * 3], vertex_normals[i * 3 + 1], vertex_normals[i * 3 + 2]);

        // std::cout << a.x << " " << a.y << " " << a.z << std::endl;
        // std::cout << b.x << " " << b.y << " " << b.z << std::endl;
        // std::cout << c.x << " " << c.y << " " << c.z << std::endl;
        // std::cout << n.x << " " << n.y << " " << n.z << std::endl << std::endl;
    }


    // Generate VAO object
    glGenVertexArrays(1, &shape_VAO);
    glBindVertexArray(shape_VAO);

    // Create VBOs for the VAO
    GLuint vertex_VBO; // Vertex Buffer
    glGenBuffers(1, &vertex_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertex_indices.size() * 3, shape_vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vVertex_attrib);
    glVertexAttribPointer(vVertex_attrib, 3, GL_FLOAT, GL_FALSE, 00, 0);
    delete[] shape_vertices;

    GLuint normal_VBO; // Normal Buffer
    glGenBuffers(1, &normal_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, normal_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertex_indices.size() * 3, vertex_normals, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vNormal_attrib);
    glVertexAttribPointer(vNormal_attrib, 3, GL_FLOAT, GL_FALSE, 0, 0);
    delete[] vertex_normals;

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}


void createVoxelizedMeshObject(unsigned int &program, unsigned int &shape_VAO, int grid_size)
{
    vector<int> vertex_indices, uv_indices, normal_indices;
    vector<glm::vec3> temp_vertices;
    vector<glm::vec2> temp_uvs;
    vector<glm::vec3> temp_normals;
    int ct = 0;

    // scale = 30; // Change Scale of the model as needed
    scale = 0.05;
    // smpl.obj takes ~ 30 as scale, buddha and bunny take ~ 0.05 as scale
    // FILE *file = fopen("src/smpl.obj", "r");
    FILE *file = fopen("src/bunny.obj", "r");
    // FILE * file = fopen("../Cobblestones3/Files/untitled.obj", "r");
    if (file == NULL)
        printf("File not found\n");

    while (true)
    {

        char head[128];
        // read the first word of the line
        int res = fscanf(file, "%s", head);
        if (res == EOF)
            break; // EOF = End Of File. Quit the loop.

        if (strcmp(head, "v") == 0)
        {
            glm::vec3 vertex;
            fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
            // printf("%f %f %f\n", vertex.x, vertex.y, vertex.z );
            vertex.x = vertex.x * scale / 10;
            vertex.y = vertex.y * scale / 10;
            vertex.z = vertex.z * scale / 10;
            temp_vertices.push_back(vertex);
            ct++;
        }

        else if (strcmp(head, "vt") == 0)
        {
            glm::vec2 uv;
            fscanf(file, "%f %f\n", &uv.x, &uv.y);
            // printf("%f %f\n", uv.x, uv.y );
            temp_uvs.push_back(uv);
        }

        else if (strcmp(head, "vn") == 0)
        {
            glm::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
            // printf("%f %f %f\n", normal.x, normal.y, normal.z );
            temp_normals.push_back(normal);
        }

        else if (strcmp(head, "f") == 0)
        {
            string vertex1, vertex2, vertex3;
            int vertexIndex[3], uvIndex[3], normalIndex[3];
            int matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2]);
            if (matches != 9)
            {
                printf("OBJ File may not contain texture coordinates or normal coordinates\n");
            }
            vertex_indices.push_back(vertexIndex[0]);
            vertex_indices.push_back(vertexIndex[1]);
            vertex_indices.push_back(vertexIndex[2]);
            uv_indices.push_back(uvIndex[0]);
            uv_indices.push_back(uvIndex[1]);
            uv_indices.push_back(uvIndex[2]);
            normal_indices.push_back(normalIndex[0]);
            normal_indices.push_back(normalIndex[1]);
            normal_indices.push_back(normalIndex[2]);
        }
    }
    fclose(file);

    glUseProgram(program);

    normalize_vertices(temp_vertices);
    // vert_stats(temp_vertices);
    voxelization(temp_vertices, grid_size);
    // create_distance_field(grid_size);
    vector<vert> vertices = create_vertices_from_voxel_grid(grid_size);
    

    // for (int i=0; i<grid_size; i++) {
    //     for (int j=0; j<grid_size; j++) {
    //         for (int k=0; k<grid_size; k++) {
    //             if (voxel_grid[i * grid_size * grid_size + j * grid_size + k] == 1) {
    //                 for (int idx=0; idx<vertices.size(); idx++) {
    //                     if (vertices[idx].pos.x == i && vertices[idx].pos.y == j && vertices[idx].pos.z == k) {
    //                         cout << vertices[idx].pos.x << " " << vertices[idx].pos.y << " " << vertices[idx].pos.z << endl;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // float min_x = 1000000, min_y = 1000000, min_z = 1000000;
    // for (int i=0; i<vertices.size(); i++) {
    //     if (vertices[i].pos.x < min_x) {
    //         min_x = vertices[i].pos.x;
    //     }
    //     if (vertices[i].pos.y < min_y) {
    //         min_y = vertices[i].pos.y;
    //     }
    //     if (vertices[i].pos.z < min_z) {
    //         min_z = vertices[i].pos.z;
    //     }
    // }
    // cout << "min_x: " << min_x << endl;
    // cout << "min_y: " << min_y << endl;
    // cout << "min_z: " << min_z << endl;

    // float max_x = -1000000, max_y = -1000000, max_z = -1000000;
    // for (int i=0; i<vertices.size(); i++) {
    //     if (vertices[i].pos.x > max_x) {
    //         max_x = vertices[i].pos.x;
    //     }
    //     if (vertices[i].pos.y > max_y) {
    //         max_y = vertices[i].pos.y;
    //     }
    //     if (vertices[i].pos.z > max_z) {
    //         max_z = vertices[i].pos.z;
    //     }
    // }

    // cout << "max_x: " << max_x << endl;
    // cout << "max_y: " << max_y << endl;
    // cout << "max_z: " << max_z << endl;

    // int voxel_count = 0;
    // for (int i=0; i<grid_size; i++) {
    //     for (int j=0; j<grid_size; j++) {
    //         for (int k=0; k<grid_size; k++) {
    //             if (voxel_grid[i * grid_size * grid_size + j * grid_size + k] == 1) {
    //                 voxel_count++;
    //             }
    //         }
    //     }
    // }

    // cout << "voxel_count: " << voxel_count << endl;

    // cout << vertices.size()/36 << endl;


    // exit(0);

    // Bind shader variables
    int vVertex_attrib = glGetAttribLocation(program, "vVertex");
    if (vVertex_attrib == -1)
    {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    int vNormal_attrib = glGetAttribLocation(program, "vNormal");
    if (vNormal_attrib == -1)
    {
        std::cout << "Could not bind location: vNormal\n";
    }


    GLfloat *shape_vertices = new GLfloat[vertices.size() * 3];
    GLfloat *vertex_normals = new GLfloat[vertices.size() * 3];

    nVertices = vertices.size() * 3;

    for (int i = 0; i < vertices.size(); i++)
    {
        shape_vertices[i * 3] = vertices[i].pos.x;
        shape_vertices[i * 3 + 1] = vertices[i].pos.y;
        shape_vertices[i * 3 + 2] = vertices[i].pos.z;
    }

    for (int i = 0; i < vertices.size(); i++)
    {
        vertex_normals[i * 3] = vertices[i].n.x;
        vertex_normals[i * 3 + 1] = vertices[i].n.y;
        vertex_normals[i * 3 + 2] = vertices[i].n.z;
    }

    // for (int i = 0; i < vertices.size(); i += 3)
    // {

    //     glm::vec3 a = glm::vec3(shape_vertices[i * 3], shape_vertices[i * 3 + 1], shape_vertices[i * 3 + 2]);
    //     glm::vec3 b = glm::vec3(shape_vertices[(i + 1) * 3], shape_vertices[(i + 1) * 3 + 1], shape_vertices[(i + 1) * 3 + 2]);
    //     glm::vec3 c = glm::vec3(shape_vertices[(i + 2) * 3], shape_vertices[(i + 2) * 3 + 1], shape_vertices[(i + 2) * 3 + 2]);
    //     glm::vec3 n = glm::vec3(vertex_normals[i * 3], vertex_normals[i * 3 + 1], vertex_normals[i * 3 + 2]);

    //     std::cout << a.x << " " << a.y << " " << a.z << std::endl;
    //     std::cout << b.x << " " << b.y << " " << b.z << std::endl;
    //     std::cout << c.x << " " << c.y << " " << c.z << std::endl;
    //     std::cout << n.x << " " << n.y << " " << n.z << std::endl << std::endl;
    // }


    // Generate VAO object
    glGenVertexArrays(1, &shape_VAO);
    glBindVertexArray(shape_VAO);

    // Create VBOs for the VAO
    GLuint vertex_VBO; // Vertex Buffer
    glGenBuffers(1, &vertex_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertices.size() * 3, shape_vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vVertex_attrib);
    glVertexAttribPointer(vVertex_attrib, 3, GL_FLOAT, GL_FALSE, 00, 0);
    delete[] shape_vertices;

    GLuint normal_VBO; // Normal Buffer
    glGenBuffers(1, &normal_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, normal_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertices.size() * 3, vertex_normals, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vNormal_attrib);
    glVertexAttribPointer(vNormal_attrib, 3, GL_FLOAT, GL_FALSE, 0, 0);
    delete[] vertex_normals;

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}






void createVoxelizedQuadMeshObject(unsigned int &program, unsigned int &shape_VAO, int grid_size)
{
    vector<int> vertex_indices, uv_indices;
    vector<glm::vec3> temp_vertices;
    vector<glm::vec2> temp_uvs;
    vector<glm::vec3> temp_normals;
    int ct = 0;

    scale = 4; // Change Scale of the model as needed
    // smpl.obj takes ~ 30 as scale, buddha and bunny take ~ 0.05 as scale
    // FILE * file = fopen("src/smpl.obj", "r");
    FILE *file = fopen("src/Cobblestones3.obj", "r");
    if (file == NULL)
        printf("File not found\n");

    while (true)
    {

        char head[128];
        // read the first word of the line
        int res = fscanf(file, "%s", head);
        if (res == EOF)
            break; // EOF = End Of File. Quit the loop.

        if (strcmp(head, "v") == 0)
        {
            glm::vec3 vertex;
            fscanf(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);
            // printf("%f %f %f\n", vertex.x, vertex.y, vertex.z );
            temp_vertices.push_back(vertex);
            ct++;
        }

        else if (strcmp(head, "vt") == 0)
        {
            glm::vec2 uv;
            fscanf(file, "%f %f\n", &uv.x, &uv.y);
            // printf("%f %f\n", uv.x, uv.y );
            temp_uvs.push_back(uv);
        }

        else if (strcmp(head, "vn") == 0)
        {
            glm::vec3 normal;
            fscanf(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);
            // printf("%f %f %f\n", normal.x, normal.y, normal.z );
            temp_normals.push_back(normal);
        }

        else if (strcmp(head, "f") == 0)
        {
            string vertex1, vertex2, vertex3;
            int vertexIndex[4], uvIndex[4];
            int matches = fscanf(file, "%d/%d %d/%d %d/%d %d/%d\n", &vertexIndex[0], &uvIndex[0], &vertexIndex[1], &uvIndex[1], &vertexIndex[2], &uvIndex[2], &vertexIndex[3], &uvIndex[3]);
            if (matches != 8)
            {
                printf("OBJ File may not contain texture coordinates\n");
            }
            vertex_indices.push_back(vertexIndex[0]);
            vertex_indices.push_back(vertexIndex[1]);
            vertex_indices.push_back(vertexIndex[2]);
            vertex_indices.push_back(vertexIndex[0]);
            vertex_indices.push_back(vertexIndex[2]);
            vertex_indices.push_back(vertexIndex[3]);
            uv_indices.push_back(uvIndex[0]);
            uv_indices.push_back(uvIndex[1]);
            uv_indices.push_back(uvIndex[2]);
            uv_indices.push_back(uvIndex[0]);
            uv_indices.push_back(uvIndex[2]);
            uv_indices.push_back(uvIndex[3]);
        }
    }
    fclose(file);

    glUseProgram(program);

    normalize_vertices(temp_vertices);
    // vert_stats(temp_vertices);
    voxelization(temp_vertices, grid_size);
    // create_distance_field(grid_size);
    vector<vert> vertices = create_vertices_from_voxel_grid(grid_size);

    // Bind shader variables
    int vVertex_attrib = glGetAttribLocation(program, "vVertex");
    if (vVertex_attrib == -1)
    {
        fprintf(stderr, "Could not bind location: vVertex\n");
        exit(0);
    }

    int vNormal_attrib = glGetAttribLocation(program, "vNormal");
    if (vNormal_attrib == -1)
    {
        std::cout << "Could not bind location: vNormal\n";
    }

    GLfloat *shape_vertices = new GLfloat[vertex_indices.size() * 3];
    GLfloat *vertex_normals = new GLfloat[vertex_indices.size() * 3];

    nVertices = vertices.size() * 3;

    for (int i = 0; i < vertices.size(); i++)
    {
        shape_vertices[i * 3] = vertices[i].pos.x;
        shape_vertices[i * 3 + 1] = vertices[i].pos.y;
        shape_vertices[i * 3 + 2] = vertices[i].pos.z;
    }

    // for (int i = 0; i < vertices.size(); i++)
    // {
    //     vertex_normals[i * 3] = vertices[i].n.x;
    //     vertex_normals[i * 3 + 1] = vertices[i].n.y;
    //     vertex_normals[i * 3 + 2] = vertices[i].n.z;
    // }

    // // Ordered vertex array to be passed to the shader to draw the mesh
    // for (int i = 0; i < vertex_indices.size(); i++)
    // {
    //     int vertexIndex = vertex_indices[i];
    //     shape_vertices[i * 3] = temp_vertices[vertexIndex - 1][0] * scale; // x-coordinate
    //     shape_vertices[i * 3 + 1] = temp_vertices[vertexIndex - 1][1] * scale; // y-coordinate
    //     shape_vertices[i * 3 + 2] = temp_vertices[vertexIndex - 1][2] * scale; // z-coordinate
    // }

    // generate normals for the triangle mesh
    for (int i = 0; i < vertex_indices.size(); i += 3)
    {
        glm::vec3 v1 = glm::vec3(shape_vertices[(i + 1) * 3] - shape_vertices[i * 3], shape_vertices[(i + 1) * 3 + 1] - shape_vertices[i * 3 + 1], shape_vertices[(i + 1) * 3 + 2] - shape_vertices[i * 3 + 2]);
        glm::vec3 v2 = glm::vec3(shape_vertices[(i + 2) * 3] - shape_vertices[(i + 1) * 3], shape_vertices[(i + 2) * 3 + 1] - shape_vertices[(i + 1) * 3 + 1], shape_vertices[(i + 2) * 3 + 2] - shape_vertices[(i + 1) * 3 + 2]);

        // using the first two vectors of a triangle, we can calculate the normal
        glm::vec3 n = cross(v1, v2);
        n = normalize(n);


        vertex_normals[i * 3] = n.x;
        vertex_normals[(i + 1) * 3] = n.x;
        vertex_normals[(i + 2) * 3] = n.x;
        vertex_normals[i * 3 + 1] = n.y;
        vertex_normals[(i + 1) * 3 + 1] = n.y;
        vertex_normals[(i + 2) * 3 + 1] = n.y;
        vertex_normals[i * 3 + 2] = n.z;
        vertex_normals[(i + 1) * 3 + 2] = n.z;
        vertex_normals[(i + 2) * 3 + 2] = n.z;
    }

    // for (int i = 0; i < vertex_indices.size(); i += 3)
    // {

    //     glm::vec3 a = glm::vec3(shape_vertices[i * 3], shape_vertices[i * 3 + 1], shape_vertices[i * 3 + 2]);
    //     glm::vec3 b = glm::vec3(shape_vertices[(i + 1) * 3], shape_vertices[(i + 1) * 3 + 1], shape_vertices[(i + 1) * 3 + 2]);
    //     glm::vec3 c = glm::vec3(shape_vertices[(i + 2) * 3], shape_vertices[(i + 2) * 3 + 1], shape_vertices[(i + 2) * 3 + 2]);
    //     glm::vec3 n = glm::vec3(vertex_normals[i * 3], vertex_normals[i * 3 + 1], vertex_normals[i * 3 + 2]);

    //     // std::cout << a.x << " " << a.y << " " << a.z << std::endl;
    //     // std::cout << b.x << " " << b.y << " " << b.z << std::endl;
    //     // std::cout << c.x << " " << c.y << " " << c.z << std::endl;
    //     // std::cout << n.x << " " << n.y << " " << n.z << std::endl << std::endl;
    // }

    // Generate VAO object
    glGenVertexArrays(1, &shape_VAO);
    glBindVertexArray(shape_VAO);

    // Create VBOs for the VAO
    GLuint vertex_VBO; // Vertex Buffer
    glGenBuffers(1, &vertex_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, vertex_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertices.size() * 3, shape_vertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vVertex_attrib);
    glVertexAttribPointer(vVertex_attrib, 3, GL_FLOAT, GL_FALSE, 00, 0);
    delete[] shape_vertices;

    GLuint normal_VBO; // Normal Buffer
    glGenBuffers(1, &normal_VBO);
    glBindBuffer(GL_ARRAY_BUFFER, normal_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * vertices.size() * 3, vertex_normals, GL_STATIC_DRAW);
    glEnableVertexAttribArray(vNormal_attrib);
    glVertexAttribPointer(vNormal_attrib, 3, GL_FLOAT, GL_FALSE, 0, 0);
    delete[] vertex_normals;

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}