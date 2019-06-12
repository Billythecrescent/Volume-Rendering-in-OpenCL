/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#define maxSteps 500
//#define tstep 0.01f


// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

int intersectBox(float4 r_o, float4 r_d, float4 boxmin, float4 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float4 invR = (float4)(1.0f,1.0f,1.0f,1.0f) / r_d;
    float4 tbot = invR * (boxmin - r_o);
    float4 ttop = invR * (boxmax - r_o);

    // re-order intersections to find smallest and largest on each axis
    float4 tmin = min(ttop, tbot);
    float4 tmax = max(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = max(max(tmin.x, tmin.y), max(tmin.x, tmin.z));
    float smallest_tmax = min(min(tmax.x, tmax.y), min(tmax.x, tmax.z));

	*tnear = largest_tmin;
	*tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}

uint findMax(uint a, uint b, uint c)
{ 
	if (a > b){ 
		if (a > c) return a;
		else return c;
	}
	else{ 
		if (b > c) return b;
		else return c;
	}
}

uint rgbaFloatToInt(float4 rgba)
{
    rgba.x = clamp(rgba.x,0.0f,1.0f);  
    rgba.y = clamp(rgba.y,0.0f,1.0f);  
    rgba.z = clamp(rgba.z,0.0f,1.0f);  
    rgba.w = clamp(rgba.w,0.0f,1.0f);  
    return ((uint)(rgba.w*255.0f)<<24) | ((uint)(rgba.z*255.0f)<<16) | ((uint)(rgba.y*255.0f)<<8) | (uint)(rgba.x*255.0f);
}

float3 calcNormal(image3d_t volume, sampler_t volumeSampler, float4 pos)
{
	float3 normal=0;
	float foot = 6.0f;
	float sample[3][3][3];
	for	(int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < 3; k++){
				sample[i][j][k] = read_imagef(volume, volumeSampler, (float4)(pos.x+foot*(i-1),pos.y+foot*(j-1),pos.z+foot*(k-1),1.0f)).x;
				//if(sample[i][j][k] > 0.01) printf("sample: %f\n",sample[i][j][k]);
			}
		}
	}
	int x = 1; int y = 1; int z = 1;

	int gradsoption = 3;
	float dx; float dy; float dz;
	//center difference algorithm
	if (gradsoption == 0){ 
		dx = sample[x+1][y][z]-sample[x-1][y][z];
		dy = sample[x][y+1][z]-sample[x-1][y-1][z];
		dz = sample[x-1][y][z+1]-sample[x-1][y][z-1];
	}
	//three-level sobel algorithm
	if (gradsoption == 1){
		dx = sample[x+1][y-1][z-1]+3*sample[x+1][y][z-1]+sample[x+1][y+1][z-1]-sample[x-1][y-1][z-1]-3*sample[x-1][y][z-1]-sample[x-1][y+1][z-1]
			+ 3*sample[x+1][y-1][z]+6*sample[x+1][y][z]+3*sample[x+1][y+1][z]-3*sample[x-1][y-1][z]-6*sample[x-1][y][z]-3*sample[x-1][y+1][z]
			+ sample[x+1][y-1][z+1]+3*sample[x+1][y][z+1]+sample[x+1][y+1][z+1]-sample[x-1][y-1][z+1]-3*sample[x-1][y][z+1]-sample[x-1][y+1][z+1];
		dy = sample[x-1][y-1][z-1]+3*sample[x][y-1][z-1]+sample[x+1][y-1][z-1]+3*sample[x-1][y][z-1]+6*sample[x][y][z-1]+3*sample[x+1][y][z-1]+sample[x-1][y+1][z-1]+3*sample[x][y+1][z-1]+sample[x+1][y+1][z-1]
			+ 0
			- sample[x-1][y-1][z+1]-3*sample[x][y-1][z+1]-sample[x+1][y-1][z+1]-3*sample[x-1][y][z+1]-6*sample[x][y][z+1]-3*sample[x+1][y][z+1]-sample[x-1][y+1][z+1]-3*sample[x][y+1][z+1]-sample[x+1][y+1][z+1];
		dz = sample[x-1][y+1][z-1]+3*sample[x][y+1][z-1]+sample[x+1][y+1][z-1]-sample[x-1][y-1][z-1]-3*sample[x][y-1][z-1]-sample[x+1][y-1][z-1]
			+ 3*sample[x-1][y+1][z]+6*sample[x][y+1][z]+3*sample[x+1][y+1][z]-3*sample[x-1][y-1][z]-6*sample[x][y-1][z]-3*sample[x+1][y-1][z]
			+ sample[x-1][y+1][z+1]+3*sample[x][y+1][z+1]+sample[x+1][y+1][z+1]-sample[x-1][y-1][z+1]-3*sample[x][y-1][z+1]-sample[x+1][y-1][z+1];
	}
	//two-level Bezier algorithm
	if (gradsoption == 2){
		dx = ((sample[x+1][y-1][z-1]-sample[x-1][y-1][z-1])+2*(sample[x+1][y-1][z]-sample[x-1][y-1][z])+(sample[x+1][y-1][z+1]-sample[x-1][y-1][z+1])
			+ 2*(sample[x+1][y][z-1]-sample[x-1][y][z-1])+4*(sample[x+1][y][z]-sample[x-1][y][z])+2*(sample[x+1][y][z+1]-sample[x-1][y][z+1])
			+ (sample[x+1][y+1][z-1]-sample[x-1][y+1][z-1])+2*(sample[x+1][y+1][z]-sample[x-1][y+1][z])+(sample[x+1][y+1][z+1]-sample[x-1][y+1][z+1]));
		dy = ((sample[x-1][y+1][z-1]-sample[x-1][y-1][z-1])+2*(sample[x-1][y+1][z]-sample[x-1][y-1][z])+(sample[x-1][y+1][z+1]-sample[x-1][y-1][z+1])
			+ 2*(sample[x][y+1][z-1]-sample[x][y-1][z-1])+4*(sample[x][y+1][z]-sample[x][y-1][z])+2*(sample[x][y+1][z+1]-sample[x][y-1][z+1])
			+ (sample[x+1][y+1][z-1]-sample[x+1][y-1][z-1])+2*(sample[x+1][y+1][z]-sample[x+1][y-1][z])+(sample[x+1][y+1][z+1]-sample[x+1][y-1][z+1]));
		dz = ((sample[x-1][y-1][z+1]-sample[x-1][y-1][z-1])+2*(sample[x-1][y][z+1]-sample[x-1][y][z-1])+(sample[x-1][y+1][z+1]-sample[x-1][y+1][z-1])
			+ 2*(sample[x][y-1][z+1]-sample[x][y-1][z-1])+4*(sample[x][y][z+1]-sample[x][y][z-1])+2*(sample[x][y+1][z+1]-sample[x][y+1][z-1])
			+ (sample[x+1][y-1][z+1]-sample[x+1][y-1][z-1])+2*(sample[x+1][y][z+1]-sample[x+1][y][z-1])+(sample[x+1][y+1][z+1]-sample[x+1][y+1][z-1]));
	}
	//Zucker-Hummel algorithm
	if (gradsoption == 3){ 
		float k1 = sqrt(3.0f)/3;
		float k2 = sqrt(2.0f)/2;
		dx = sample[x+1][y][z]-sample[x-1][y][z]+
			k1*(sample[x+1][y+1][z+1]-sample[x-1][y+1][z+1]+sample[x+1][y-1][z+1]-sample[x-1][y-1][z+1]+
			sample[x+1][y+1][z-1]-sample[x-1][y+1][z-1]+sample[x+1][y-1][z-1]-sample[x-1][y-1][z-1])+
			k2*(sample[x+1][y][z+1]-sample[x-1][y][z+1]+sample[x+1][y+1][z+1]-sample[x-1][y+1][z+1]+
			sample[x+1][y][z-1]-sample[x-1][y][z-1]+sample[x+1][y-1][z]-sample[x-1][y-1][z]);
		dy = sample[x][y+1][z]-sample[x][y-1][z]+
			k1*(sample[x-1][y+1][z-1]-sample[x-1][y-1][z+1]+sample[x+1][y-1][z+1]-sample[x+1][y-1][z+1])+
			sample[x-1][y+1][z-1]-sample[x-1][y-1][z-1]+sample[x+1][y+1][z-1]-sample[x+1][y-1][z-1]+
			k2*(sample[x][y+1][x+1]-sample[x-1][y][z+1]+sample[x+1][y+1][z]-sample[x+1][y+1][z]+
			sample[x][y+1][z-1]-sample[x][y-1][z-1]+sample[z+1][y-1][z]-sample[x+1][y-1][z]);
		dz = sample[x][y][z+1]-sample[x][y][z-1]+
			k1*(sample[x-1][y+1][z+1]-sample[x-1][y+1][z-1]+sample[x+1][y+1][z+1]-sample[x+1][y+1][z-1]+
			sample[x-1][y-1][z+1]-sample[x-1][y-1][z-1]+sample[x+1][y-1][z+1]-sample[x+1][y-1][z-1])+
			k2*(sample[x][y+1][z+1]-sample[x][y+1][z+1]+sample[x][y-1][z+1]-sample[x][y-1][z-1]+
			sample[x+1][y][z+1]-sample[x+1][y][z-1]+sample[x-1][y][z+1]-sample[x-1][y][z-1]);
	}
	normal = normalize((float3)(dx,dy,dz));
	
	//if (normal.x * normal.y * normal.z != 0)printf("normal: [%f,%f,%f],[%f, %f, %f]\n",dx,dy,dz,normal.x,normal.y,normal.z);
	return normal;
}

float3 calcLight(float4 eyeRay_d, image3d_t volume, sampler_t volumeSampler, float4 pos)
{
	float3 light_ambient = (float3)(0.2f, 0.2f, 0.2f);
	float3 light_diffuse = (float3)(0.6f, 0.6f, 0.6f);
	float3 light_specular = (float3)(0.6f, 0.6f, 0.6f);
	float3 material_ambient = (float3)(0.9f, 0.7f, 0.5f);
	float3 material_diffuse = (float3)(0.9f, 0.7f, 0.5f);
	float3 material_specular = (float3)(0.6f, 0.6f, 0.6f);
	float shininess = 32;
	float3 V = eyeRay_d.xyz;
	float3 L = normalize((float3)(1.0f,1.0f,1.0f)); 
	float3 H = normalize(V + L); //∞Î≥ÃœÚ¡ø
	float3 N = calcNormal(volume, volumeSampler, pos);
	float3 ambientcolor = light_ambient * material_ambient;
	float3 diffusecolor = light_diffuse * material_diffuse * max(dot(N, L), 0.0f);
	float3 specularcolor = light_specular * material_specular * pow(max(dot(N, H), 0.0f), shininess);
	//float3 color = ambientcolor;
	float3 color = ambientcolor + diffusecolor + specularcolor;
	//if (dot(N,L) != 0) printf("N,L,dot(N,L): [%f,%f,%f],[%f,%f,%f], %f\n",N.x,N.y,N.z,L.x,L.y,L.z,dot(N,L));
	return color;
}

__kernel void
d_render(__global uint *d_output, 
         uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale,
         __constant float* invViewMatrix
 #ifdef IMAGE_SUPPORT
          ,__read_only image3d_t volume,
          __read_only image2d_t transferFunc,
          sampler_t volumeSampler,
          sampler_t transferFuncSampler
 #endif
         )

{	
    uint x = get_global_id(0);
    uint y = get_global_id(1);

    float u = (x / (float) imageW)*2.0f-1.0f;
    float v = (y / (float) imageH)*2.0f-1.0f;

	
	int width = get_image_width(volume);
	int height = get_image_height(volume);
	int lengt = get_image_depth(volume);
	int maxLen = findMax(lengt, width, height);
	float wid = (float)width / (float)maxLen;
	float hei = (float)height / (float)maxLen;
	float len = (float)lengt / (float)maxLen;

    //float tstep = 0.01f;

    float4 boxMin = (float4)(-wid, -hei, -len, 1.0f);
    float4 boxMax = (float4)(wid, hei, len, 1.0f);
	//float4 boxMin = (float4)(-1.0f, -1.0f, -1.0f,1.0f);
    //float4 boxMax = (float4)(1.0f, 1.0f, 1.0f,1.0f);

    // calculate eye ray in world space
    float4 eyeRay_o;
    float4 eyeRay_d;

    eyeRay_o = (float4)(invViewMatrix[3], invViewMatrix[7], invViewMatrix[11], 1.0f);   
	//printf("eyeRay_o: %f, %f, %f, %f\n",eyeRay_o.x,eyeRay_o.y,eyeRay_o.z,eyeRay_o.w);
    float4 temp = normalize(((float4)(u, v, -2.0f,0.0f)));
    eyeRay_d.x = dot(temp, ((float4)(invViewMatrix[0],invViewMatrix[1],invViewMatrix[2],invViewMatrix[3])));
    eyeRay_d.y = dot(temp, ((float4)(invViewMatrix[4],invViewMatrix[5],invViewMatrix[6],invViewMatrix[7])));
    eyeRay_d.z = dot(temp, ((float4)(invViewMatrix[8],invViewMatrix[9],invViewMatrix[10],invViewMatrix[11])));
    eyeRay_d.w = 0.0f;

    // find intersection with box
	float tnear, tfar;
	int hit = intersectBox(eyeRay_o, eyeRay_d, boxMin, boxMax, &tnear, &tfar);
    if (!hit) {
        if ((x < imageW) && (y < imageH)) {
            // write output color
            uint i =(y * imageW) + x;
            d_output[i] = 0;
        }
        return;
    }
	if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from back to front, accumulating color
    temp = (float4)(0.0f,0.0f,0.0f,0.0f);
	float a = 0.0f;
	float tstep = 0.01f;
    float t = tfar;
	float accumA = 0.0f;
	float3 accumC = (float3)(0.0f,0.0f,0.0f);
    for(uint i=0; i<maxSteps; i++) {		
        float4 pos = eyeRay_o + eyeRay_d*t;
        pos.x = (pos.x + wid) * 0.5f / wid * (width - 1);   // pos.x ranged in (-wid, wid)
		pos.y = (pos.y + hei) * 0.5f / hei * (height - 1);
		pos.z = (pos.z + len) * 0.5f / len * (lengt - 1);
		//pos = pos*0.5f+0.5f;

        // read from 3D texture        
#ifdef IMAGE_SUPPORT        
        float4 sample = read_imagef(volume, volumeSampler, pos);
        //printf("pos: %f, %f, %f, %f\n", pos.x,pos.y,pos.z,pos.w);
		
		// lookup in transfer function texture
        float2 transfer_pos = (float2)((sample.x-transferOffset)*transferScale, 0.5f);
		//printf("transpos: %f, %f\n", transfer_pos.x,transfer_pos.y);
        float4 col = read_imagef(transferFunc, transferFuncSampler, transfer_pos);
#else
        float4 col = (float4)(pos.x,pos.y,pos.z,.25f);
#endif

        // accumulate result
		//col = (float4)calcLight(eyeRay_d, volume, volumeSampler, pos);
		float a = col.w*density;
		//if (col.x > 0.2) printf("col: %f, %f, %f, %f\n",col.x,col.y,col.z,col.w);
		float3 color = (float3)calcLight(eyeRay_d, volume, volumeSampler, pos);
		//float3 color = (float3)(col.x,col.y,col.z);
		accumC = a * color + (1-a)*accumC;		
		accumA = (1 - accumA) * a + accumA;

		col.xyz = accumC;
		col.w = accumA;
        temp = mix(temp, col, (float4)(a, a, a, a));
		//printf("R = %d\n",temp[0]);
		//if (temp.x*temp.y*temp.z != 0) printf("temp: %f, %f, %f, %f\n",temp.x,temp.y,temp.z,temp.w);

		tstep = 0.00001*sqrt(pow(pos.x-eyeRay_d.x,2) + pow(pos.y-eyeRay_d.y,2) + pow(pos.z-eyeRay_d.z,2));
        t -= tstep;
        if (t < tnear) break;
    }
    temp *= brightness;

    if ((x < imageW) && (y < imageH)) {
        // write output color
        uint i =(y * imageW) + x;
        d_output[i] = rgbaFloatToInt(temp);
    }
	
}

