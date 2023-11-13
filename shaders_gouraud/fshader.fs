#version 330 core

in vec3 n;
in vec3 e;
in vec3 l;
out vec4 outColor;

vec3 Ls = vec3(1.0, 1.0, 1.0);
vec3 Ld = vec3(0.7, 0.7, 0.7);
vec3 La = vec3(0.6, 0.3, 0.4);

vec3 ks = vec3(1.0, 1.0, 1.0);
vec3 kd = vec3(0.5, 0.6, 0.4);
vec3 ka = vec3(1.0, 1.0, 1.0);
// vec3 ka = vec3(0.0, 0.0, 0.0);

float spec_exp = 5.0;
// with h, it looks like spec_exp is 1.0

//ambient
// vec3 Ia = ka * La;
vec3 Ia = vec3(0.0, 0.0, 0.0);

//diffuse
vec3 Id = kd * max(dot(n, l) * Ld, 0.0);
// vec3 Id = vec3(0.0, 0.0, 0.0);

//specular
vec3 r = 2.0 * dot(n, l) * n - l;
// e = normalize(e);
vec3 h = normalize(l - e);

// - to keep both diffuse and specular light direction same
// vec3 Is = ks * pow (max( dot(n, h), 0.0), spec_exp) * Ls;
// vec3 Is = ks * pow( max( dot(r, -e), 0.0), spec_exp) * Ls;
vec3 Is = vec3(0.0, 0.0, 0.0);


vec3 fColor = Ia + Id +Is;
// vec3 fColor = vec3(0.0, 0.0, 0.0);

void main(void) {
        outColor = vec4(fColor, 1.0);
}
