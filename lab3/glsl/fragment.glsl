varying vec3 normal;
varying vec4 world_position;

void main()
{
    vec4 ambient = vec4(1, 0, 0, 1);
    vec4 diffuse = vec4(0, 1, 0, 1);
    vec4 specular = vec4(0, 0, 1, 1);

    vec3 l = vec3(gl_LightSource[0].position.xyz - world_position.xyz);
    vec3 L = normalize(l);
    vec3 n = normalize(normal);

    ambient = gl_LightModel.ambient * gl_FrontMaterial.ambient;

    float diffuse_intensity = max(dot(L, n),0.0);
    diffuse = gl_FrontMaterial.diffuse * gl_LightSource[0].diffuse * diffuse_intensity;


    vec3 r = 2.0 * dot(L, normal) * normal - L;
    float specular_intensity = pow(max(dot(r, -normalize(world_position.xyz)),0.0),gl_FrontMaterial.shininess);
    specular = gl_FrontMaterial.specular * specular_intensity;
    
    gl_FragColor = ambient + diffuse + specular;
}
