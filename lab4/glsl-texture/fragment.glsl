uniform sampler2D tex;
varying vec3 N;
varying vec4 position;

void main()
{
    //vec3 diff = vec3(gl_LightSource[0].position - position);
    //vec3 L = normalize(diff);
    //vec3 normal = normalize(N);

    //float intensity = max(dot(L, normal),0.0);
    //gl_FragColor.rgb = intensity * gl_FrontMaterial.diffuse.rgb * gl_LightSource[0].diffuse.rgb + gl_FrontMaterial.ambient.rgb * //gl_LightSource[0].ambient.rgb;
    //gl_FragColor.a = 1.0;


    // get the texture color from "tex", using a texture lookup function 
    // and the s,t coordinates of
    vec4 tex_color = (texture2D(tex,gl_TexCoord[0].st));

    // set gl_FragColor to tex_color
    gl_FragColor = tex_color;

}
