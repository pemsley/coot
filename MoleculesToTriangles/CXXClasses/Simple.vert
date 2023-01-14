uniform float Scale;

void main(void)
{
    vec4 a = gl_Vertex;
    a.x = a.x * Scale;
    a.y = a.y * Scale;
    gl_Position = gl_ModelViewProjectionMatrix * a;
    
}
