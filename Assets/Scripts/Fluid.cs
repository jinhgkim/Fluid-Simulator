using UnityEngine;

public class Fluid : MonoBehaviour
{
    const int size = 256;

    float dt; // timestep
    float diff; // diffusion
    float visc; // viscosity

    float[] s;
    float[] density;

    float[] vx;
    float[] vy;

    float[] vx0;
    float[] vy0;

    static int IX(int x, int y) => x + y * size;

    public Fluid(float dt, float diff, float visc)
    {
        this.dt = dt;
        this.diff = diff;
        this.visc = visc;

        int area = size * size;

        this.s = new float[area];
        this.density = new float[area];

        this.vx = new float[area];
        this.vy = new float[area];

        this.vx0 = new float[area];
        this.vy0 = new float[area];
    }

    void AddDensity(int x, int y, float amount)
    {
        int index = IX(x, y);

        density[index] += amount;
    }

    void AddVelocity(int x, int y, float amountX, float amountY)
    {
        int index = IX(x, y);

        vx[index] += amountX;
        vy[index] += amountY;
    }


}
