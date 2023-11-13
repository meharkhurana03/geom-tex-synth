#version 330 core

in vec3 n;
in vec3 ePos;
in vec3 lPos;
in vec3 fPos;
in vec3 fColor;

out vec4 outColor;

vec3 Ls = vec3(1.0, 1.0, 1.0);
vec3 Ld = vec3(0.7, 0.7, 0.7);
vec3 La = vec3(0.7, 0.3, 0.3);

vec3 ks = vec3(1.0, 1.0, 1.0);
vec3 kd = vec3(0.5, 0.6, 0.4);
vec3 ka = vec3(1.0, 1.0, 1.0);
// vec3 ka = vec3(0.0, 0.0, 0.0);

float spec_exp = 5.0;

//ambient
vec3 Ia = ka * La;
// vec3 Ia = vec3(0.0, 0.0, 0.0);

vec3 l = normalize(lPos - fPos);
vec3 e = normalize(ePos - fPos);
vec3 n_norm = normalize(n);
vec3 r = normalize(reflect(l, n_norm));

//diffuse
vec3 Id = kd * max(dot(n_norm, l) * Ld, 0.0);
// vec3 Id = vec3(0.0, 0.0, 0.0);

//specular
// vec3 Is = ks * pow ( max(dot(r, e), 0.0), spec_exp ) * Ls;

// - to keep both diffuse and specular light direction same
vec3 h = normalize(l - e);
vec3 Is = ks * max( pow(dot(n, h), spec_exp), 0.0) * Ls;
// vec3 Is = vec3(0.0, 0.0, 0.0);


vec3 color = (Ia + Id + Is) * fColor;
// vec3 fColor = vec3(0.0, 0.0, 0.0);

void main(void) {
        outColor = vec4(color, 1.0);
}


// #version 330 core

// in vec3 n;
// in vec3 ePos;
// in vec3 lPos;
// in vec3 fPos;
// in vec3 fColor;
// out vec4 outColor;

// vec3 Ls = vec3(0.9, 0.9, 0.9);
// vec3 Ld = vec3(0.8, 0.8, 0.8);
// vec3 La = vec3(0.7, 0.3, 0.3);

// vec3 ks = vec3(0.7, 0.7, 0.7);
// vec3 kd = vec3(0.5, 0.6, 0.4);
// vec3 ka = vec3(1.0, 1.0, 1.0);

// float spec_exp = 5.0;

// vec3 norm_l = normalize(lPos - fPos);

// //ambient
// vec3 Ia = ka * La;

// vec3 norm_n = normalize(n);

// //diffuse
// vec3 Id = kd * max(dot(norm_n, norm_l) * Ld, 0.0);

// //specular
// vec3 norm_e = normalize(ePos - fPos);
// // vec3 norm_e = normalize(ePos);

// // - to keep both diffuse and specular light direction same
// // vec3 h = normalize(norm_l - norm_e);
// // vec3 Is = ks * max( pow(dot(norm_n, h), spec_exp), 0.0) * Ls;

// vec3 r = normalize(reflect(norm_l, norm_n));
// vec3 Is = ks * pow ( max(dot(r, norm_e), 0.0), spec_exp ) * Ls;


// vec3 color = (Ia + Id + Is) * fColor;

// void main(void) {
//         outColor = vec4(color, 1.0);
// }
