/************************************************************************
     File:        MazeWindow.cpp

     Author:     
                  Stephen Chenney, schenney@cs.wisc.edu
     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu

     Comment:    
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for the MazeWindow class. The MazeWindow is
						the window in which the viewer's view of the maze is displayed.
		

     Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "MazeWindow.h"
#include <Fl/math.h>
#include <Fl/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <glm/glm.hpp>


//*************************************************************************
//
// * Constructor
//=========================================================================
MazeWindow::
MazeWindow(int x, int y, int width, int height, const char *label,Maze *m)
	: Fl_Gl_Window(x, y, width, height, label)
//=========================================================================
{
	maze = m;

	// The mouse button isn't down and there is no key pressed.
	down = false;
	z_key = 0;
}


//*************************************************************************
//
// * Set the maze. Also causes a redraw.
//=========================================================================
void MazeWindow::
Set_Maze(Maze *m)
//=========================================================================
{
	// Change the maze
	maze = m;

	// Force a redraw
	redraw();
}


//*************************************************************************
//
// * draw() method invoked whenever the view changes or the window
//   otherwise needs to be redrawn.
//=========================================================================
void MazeWindow::
draw(void)
//=========================================================================
{
	float   focal_length;

	if ( ! valid() ) {
		// The OpenGL context may have been changed
		// Set up the viewport to fill the window.
		glViewport(0, 0, w(), h());

		// We are using orthogonal viewing for 2D. This puts 0,0 in the
		// middle of the screen, and makes the image size in view space
		// the same size as the window.
		gluOrtho2D(-w() * 0.5, w() * 0.5, -h() * 0.5, h() * 0.5);

		// Sets the clear color to black.
		glClearColor(0.0, 0.0, 0.0, 1.0);
	}

	// Clear the screen.
	glClear(GL_COLOR_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glBegin(GL_QUADS);
		// Draw the "floor". It is an infinite plane perpendicular to
		// vertical, so we know it projects to cover the entire bottom
		// half of the screen. Walls of the maze will be drawn over the top
		// of it.
		glColor3f(0.2f, 0.2f, 0.2f);
		glVertex2f(-w() * 0.5f, -h() * 0.5f);
		glVertex2f( w() * 0.5f, -h() * 0.5f);
		glVertex2f( w() * 0.5f, 0.0       );
		glVertex2f(-w() * 0.5f, 0.0       );

		// Draw the ceiling. It will project to the entire top half
		// of the window.
		glColor3f(0.4f, 0.4f, 0.4f);
		glVertex2f( w() * 0.5f,  h() * 0.5f);
		glVertex2f(-w() * 0.5f,  h() * 0.5f);
		glVertex2f(-w() * 0.5f, 0.0       );
		glVertex2f( w() * 0.5f, 0.0       );
	glEnd();


	if ( maze ) {
		// Set the focal length. We can do this because we know the
		// field of view and the size of the image in view space. Note
		// the static member function of the Maze class for converting
		// radians to degrees. There is also one defined for going backwards.
		focal_length = w() / (float)(2.0*tan(Maze::To_Radians(maze->viewer_fov)*0.5));
	
		// Draw the 3D view of the maze (the visible walls.) You write this.
		// Note that all the information that is required to do the
		// transformations and projection is contained in the Maze class,
		// plus the focal length.

		glClear(GL_DEPTH_BUFFER_BIT);
		
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		float aspect = (float)w() / h();
		float f = tan(glm::radians(maze->viewer_fov) / 2.0f);
		f = 1.0f / f;
		float zNear = 0.01f;
		float zFar = 200.0f;
		float perspectiveMatrix[16] = {
			f / aspect, 0, 0, 0,
			0, f, 0, 0, 
			0, 0, (zFar + zNear) / (zNear - zFar), -1,
			0, 0, (2 * zFar * zNear) / (zNear - zFar), 0
		};
		//gluPerspective(maze->viewer_fov, aspect, 0.01, 200);
		//glMultMatrixf(perspectiveMatrix);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		float viewer_pos[3] = { maze->viewer_posn[Maze::Y], 0.0f, maze->viewer_posn[Maze::X] };
		float eye[3] = { viewer_pos[Maze::X], viewer_pos[Maze::Y], viewer_pos[Maze::Z] };
		float at[3] = { viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)),
			viewer_pos[Maze::Y],
			viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)) };

		glm::vec3 forward = { sin(Maze::To_Radians(maze->viewer_dir)) , 0, cos(Maze::To_Radians(maze->viewer_dir)) };
		forward = glm::normalize(forward);
		glm::vec3 up = { 0, 1, 0 };
		glm::vec3 right = glm::cross(forward, up);
		up = glm::cross(right, forward);

		float lookat[] = {
			right.x, up.x, -forward.x, 0,
			right.y, up.y, -forward.y, 0,
			right.z, up.z, -forward.z, 0,
			-(right.x * viewer_pos[Maze::X] + right.y * viewer_pos[Maze::Y] + right.z * viewer_pos[Maze::Z]),
			-(up.x * viewer_pos[Maze::X] + up.y * viewer_pos[Maze::Y] + up.z * viewer_pos[Maze::Z]),
			(forward.x * viewer_pos[Maze::X] + forward.y * viewer_pos[Maze::Y] + forward.z * viewer_pos[Maze::Z]),
			1.0f
		};

		/*
		gluLookAt(viewer_pos[Maze::X], viewer_pos[Maze::Y], viewer_pos[Maze::Z],
			viewer_pos[Maze::X] + sin(Maze::To_Radians(maze->viewer_dir)),
			viewer_pos[Maze::Y],
			viewer_pos[Maze::Z] + cos(Maze::To_Radians(maze->viewer_dir)),
			0.0, 1.0, 0.0);*/
		
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		//glMultMatrixf(lookat);

		maze->Draw_View(focal_length, lookat, perspectiveMatrix);
	}
}


//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Drag(float dt)
//=========================================================================
{
	float   x_move, y_move, z_move;

	if ( down ) {
		int	dx = x_down - x_last;
		int   dy = y_down - y_last;
		float dist;

		// Set the viewing direction based on horizontal mouse motion.
		maze->Set_View_Dir(d_down + 360.0f * dx / (float)w());

		// Set the viewer's linear motion based on a speed (derived from
		// vertical mouse motion), the elapsed time and the viewing direction.
		dist = 10.0f * dt * dy / (float)h();
		x_move = dist * (float)cos(Maze::To_Radians(maze->viewer_dir));
		y_move = dist * (float)sin(Maze::To_Radians(maze->viewer_dir));
	}
	else {
		x_move = 0.0;
		y_move = 0.0;
	}

	// Update the z posn
	z_move = z_key * 0.1f;
	z_key = 0;

	// Tell the maze how much the view has moved. It may restrict the motion
	// if it tries to go through walls.
	maze->Move_View_Posn(x_move, y_move, z_move);

	return true;
}

//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
KeyEvent(float dt)
//=========================================================================
{
	float speed = 0.1f;
	maze->Set_View_Dir(maze->viewer_dir + d_keyDown);
	maze->Move_View_Posn(speed * x_dirMove, speed * y_dirMove, 0);
	
	return true;
}
//*************************************************************************
//
// *
//=========================================================================
bool MazeWindow::
Update(float dt)
//=========================================================================
{
	// Update the view

	if ( down || z_key ) // Only do anything if the mouse button is down.
		return Drag(dt);
	if ( keyDown )
		//return KeyEvent(dt);
	
	// Nothing changed, so no need for a redraw.
	return false;
}


//*************************************************************************
//
// *
//=========================================================================
int MazeWindow::
handle(int event)
//=========================================================================
{
	if (!maze)
		return Fl_Gl_Window::handle(event);

	// Event handling routine.
	switch ( event ) {
		case FL_PUSH:
			down = true;
			x_last = x_down = Fl::event_x();
			y_last = y_down = Fl::event_y();
			d_down = maze->viewer_dir;
			return 1;
		case FL_DRAG:
			x_last = Fl::event_x();
			y_last = Fl::event_y();
			return 1;
			case FL_RELEASE:
			down = false;
			return 1;
		case FL_KEYDOWN:
			keyDown = true;
			if (Fl::event_key() == 'w') {
				x_dirMove = (float)cos(Maze::To_Radians(maze->viewer_dir));
				y_dirMove = (float)sin(Maze::To_Radians(maze->viewer_dir));
				return 1;
			}
			if (Fl::event_key() == 's') {
				x_dirMove = -(float)cos(Maze::To_Radians(maze->viewer_dir));
				y_dirMove = -(float)sin(Maze::To_Radians(maze->viewer_dir));
				return 1;
			}
			if (Fl::event_key() == 'a') {
				x_dirMove = (float)cos(Maze::To_Radians((double)maze->viewer_dir + 90.0));
				y_dirMove = (float)sin(Maze::To_Radians((double)maze->viewer_dir + 90.0));
				return 1;
			}
			if (Fl::event_key() == 'd') {
				x_dirMove = (float)cos(Maze::To_Radians((double)maze->viewer_dir - 90.0));
				y_dirMove = (float)sin(Maze::To_Radians((double)maze->viewer_dir - 90.0));
				return 1;
			}
			if (Fl::event_key() == FL_Right) {
				d_keyDown = -2.0f;
				return 1;
			}
			if (Fl::event_key() == FL_Left) {
				d_keyDown = 2.0f;
				return 1;
			}
			return Fl_Gl_Window::handle(event);
		case FL_KEYUP:
			keyDown = false;
			x_dirMove = 0.0f;
			y_dirMove = 0.0f;
			d_keyDown = 0.0f;
			return 1;
		case FL_FOCUS:
		case FL_UNFOCUS:
			return 1;
	}

	// Pass any other event types on the superclass.
	return Fl_Gl_Window::handle(event);
}


