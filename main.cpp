#include <deque>
#include <raylib.h>
#include <cmath>

constexpr unsigned int window_width{1200};
constexpr unsigned int window_height{900};
const double cx = static_cast<double>(window_width) / 2.0;
const double cy = static_cast<double>(window_height) / 3.0;
constexpr double g{981};

  
  
struct State{
    double alpha1, beta1;
    double alpha2, beta2;

};

  

State computederivative(const State& s, double m1, double m2, double l1, double l2){
    State derived;
    derived.alpha1 = s.beta1;
    derived.alpha2 = s.beta2;

    double delta{s.alpha1-s.alpha2};
    double denom = (2*m1 + m2 - m2* cos(2*s.alpha1-2*s.alpha2));

    derived.beta1 = (-g *(2*m1 + m2)*sin(s.alpha1) - m2*g*sin(s.alpha1-s.alpha2*2)- 2*sin(delta)*m2*(l2*s.beta2*s.beta2 + l1*s.beta1*s.beta1 * cos(delta))) / (l1*denom);
    derived.beta2 = (2*sin(delta)*(l1*s.beta1*s.beta1*(m1+m2) + g*(m1+m2) * cos(s.alpha1)+ s.beta2*s.beta2*l2*m2*cos(delta))) / (l2*denom);

    return derived;

}


void rk4step(State& s, double dt, double m1, double m2, double l1, double l2){

    State k1{computederivative(s, m1,m2,l1,l2)};

    State s2{
        s.alpha1 + 0.5*dt * k1.alpha1,
        s.beta1 + 0.5*dt * k1.beta1,
        s.alpha2 + 0.5*dt * k1.alpha2,
        s.beta2 + 0.5*dt * k1.beta2,
    };

    State k2{computederivative(s2, m1,m2,l1,l2)};

    State s3{
        s.alpha1 + 0.5*dt * k2.alpha1,
        s.beta1 + 0.5*dt * k2.beta1,
        s.alpha2 + 0.5*dt * k2.alpha2,
        s.beta2 + 0.5*dt * k2.beta2,
    };

    State k3{computederivative(s3, m1,m2,l1,l2)};

    State s4{
        s.alpha1 + dt * k3.alpha1,
        s.beta1 + dt * k3.beta1,
        s.alpha2 + dt * k3.alpha2,
        s.beta2 + dt * k3.beta2,
    };

    State k4{computederivative(s4, m1,m2,l1,l2)};

    s.alpha1 += (dt/6.0) * (k1.alpha1 + k2.alpha1*2 + k3.alpha1*2 + k4.alpha1);
    s.beta1 += (dt/6.0) * (k1.beta1 + k2.beta1*2 + k3.beta1*2 + k4.beta1);
    s.alpha2 += (dt/6.0) * (k1.alpha2 + k2.alpha2*2 + k3.alpha2*2 + k4.alpha2);
    s.beta2 += (dt/6.0) * (k1.beta2 + k2.beta2*2 + k3.beta2*2 + k4.beta2);
    
}

int main(){
    InitWindow(window_width, window_height, "pend");
    SetTargetFPS(60);

    constexpr double m1{30};
    constexpr double m2{20};
    constexpr double l1{200};
    constexpr double l2{200};
    constexpr double dt{0.02};

    constexpr int maxTrailLength = 150;
    std::deque<Vector2> trail{};

    State startstate{PI/2, 0,PI/4, 0};

    while(!WindowShouldClose()){
        rk4step(startstate, dt, m1, m2,l1,l2);

        double alpha1_rad = startstate.alpha1;
        double alpha2_rad = startstate.alpha2;

        double x1 = cx + l1 * sin(alpha1_rad);
        double y1 = cy + l1 * cos(alpha1_rad);
        double x2 = x1 + l2 * sin(alpha2_rad);
        double y2 = y1 + l2 * cos(alpha2_rad);

        trail.push_back(Vector2{static_cast<float>(x2),static_cast<float>(y2)});
        while (trail.size() > maxTrailLength) {trail.pop_front();}

        for(int i{1}; i<trail.size(); ++i){
            DrawLineEx(trail[i-1], trail[i], 3, Fade(SKYBLUE, static_cast<float>(i) / trail.size()));
        }

    

        BeginDrawing();
        ClearBackground(BLACK);

        DrawLine(cx, cy, x1, y1, WHITE);
        DrawLine(x1, y1, x2, y2, WHITE);
        
        
        

        DrawCircle(static_cast<int>(x1), static_cast<int>(y1), 14.0f, RED);
        DrawCircle(static_cast<int>(x2), static_cast<int>(y2), 14.0f, BLUE);
       

        EndDrawing();
    }
    CloseWindow();
    return 0;
}
