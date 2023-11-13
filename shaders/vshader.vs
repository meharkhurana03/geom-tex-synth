#version 330 core

in vec3 vVertex;
in vec3 vNormal;

uniform mat4 vModel;
uniform mat4 vView;
uniform mat4 vProjection;
uniform vec3 lpos_world;
uniform vec3 eye_normal;
uniform vec3 vColor;

out vec3 n;
out vec3 ePos;
out vec3 lPos;
out vec3 fColor;
out vec3 fPos;

void main() {
	gl_Position = vProjection * vView * vModel * vec4(vVertex, 1.0);
    fPos = vec3(vModel * vec4(vVertex, 1.0));
    n = vNormal;
	// n = normalize(vNormal);
    lPos = lpos_world;
    ePos = eye_normal;
    // e = normalize(eye_normal);
    fColor = vColor;
}


// #version 330 core

// in vec3 vVertex;
// in vec3 vNormal;

// uniform mat4 vModel;
// uniform mat4 vView;
// uniform mat4 vProjection;
// uniform vec3 vColor;
// uniform vec3 lpos_world;
// uniform vec3 eye_normal;


// out vec3 n;
// out vec3 ePos;
// out vec3 lPos;
// out vec3 fColor;
// out vec3 fPos;


// void main() {
// 	gl_Position = vProjection * vView * vModel * vec4(vVertex, 1.0);
// 	fPos = vec3(vModel * vec4(vVertex, 1.0));
// 	fColor = vColor; //Interpolate color
//     n = vNormal;
	
// 	lPos = lpos_world;
//     ePos = eye_normal;
	
// }