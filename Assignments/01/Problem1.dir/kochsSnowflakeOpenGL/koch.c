#include <GL/glut.h>
#include <math.h>

GLfloat oldx=-0.7,oldy=0.5;
GLint myIterations = 6;
GLfloat myLength ( GLint iter ) {
    return 0.05 * pow(3,3) / pow ( 3, iter );
}

void drawkoch(GLfloat dir,GLfloat len,GLint iter) {
    GLdouble dirRad = 0.0174533 * dir;  
    GLfloat newX = oldx + len * cos(dirRad);
    GLfloat newY = oldy + len * sin(dirRad);
    if (iter==0) {
        glVertex2f(oldx, oldy);
        glVertex2f(newX, newY);
        oldx = newX;
        oldy = newY;
    }
    else {
        iter--;
        //draw the four parts of the side _/\_
        drawkoch(dir, len, iter);
        dir += 60.0;
        drawkoch(dir, len, iter);
        dir -= 120.0;
        drawkoch(dir, len, iter);
        dir += 60.0;
        drawkoch(dir, len, iter);
    }
}

void mydisplay(){
    glClear( GL_COLOR_BUFFER_BIT );
    glBegin(GL_LINES);
    glColor3f(0.0, 0.54, 1.0); // make it red

    drawkoch(0.0,myLength(myIterations),myIterations);
    drawkoch(-120.0,myLength(myIterations), myIterations);
    drawkoch(120.0,myLength(myIterations),myIterations);

    glEnd();
    glFlush();
}

int main(int argc, char** argv)
{
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);
    glutInitWindowSize(500,500);
    glutInitWindowPosition(0,0);
    glutCreateWindow("The ting goes SKRRRAH! Pr pim pr pim pim. BOOM.");
    glutDisplayFunc(mydisplay);
    glutMainLoop();
}
