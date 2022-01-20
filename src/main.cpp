#include "Fluid.h"
#include <SFML/Graphics.hpp>
#include <SFML/Graphics/Image.hpp>

#include "DummyFluid.h"

int main()
{
	//fluid parameters
	int xSize, ySize ;
	xSize = ySize = 50;
	int elementSize = 20;
	float diff = 0.107;
	float visc = 0.19975;
	float dens_value = 100;
	//simulation parameters
	int framerate = 60;
	float dt = 1.0 / framerate;
	int numIterations = 50;
	float eps = 1;
	//mouse parameters
	int mouseX = 0;
	int mouseY = 0;
	int lastMouseX = 0;
	int lastMouseY = 0;

	//setting up render window
	sf::RenderWindow window(sf::VideoMode((xSize) * elementSize + 1, (ySize) * elementSize + 1), "Project");
	window.setFramerateLimit(framerate);

	//creating rendering tools
	sf::Image image;
	sf::Texture texture;
	sf::Sprite sprite;
	image.create((xSize) * elementSize + 1, (ySize) * elementSize + 1, sf::Color::Black);

	//now let's create the fluid we want (x, y, diffusion, viscosity, timestep, accuracy/quality)
	 auto fluid = new FluidCube(ySize, diff, visc, dt, numIterations);
	while (window.isOpen())
	{
		//window details
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
				window.close();
		}

		lastMouseX = mouseX;
		lastMouseY = mouseY;
		sf::Vector2i localPosition = sf::Mouse::getPosition(window);
		mouseX = localPosition.x;
		mouseY = localPosition.y;

		if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
			if ((mouseX / elementSize > 0 && mouseX / elementSize < xSize - 1 && mouseY / elementSize  > 0 && mouseY / elementSize < ySize - 1) && (lastMouseX / elementSize > 0 && lastMouseX / elementSize < xSize - 1 && lastMouseY / elementSize  > 0 && lastMouseY / elementSize < ySize - 1)) {
						 fluid->addDensity(lastMouseX / elementSize, lastMouseY / elementSize, dens_value);
			}
		}

		if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
			fluid->reset();
		}
		fluid->update();
		fluid->render(elementSize, image);
		texture.loadFromImage(image);
		sprite.setTexture(texture, true);
		window.draw(sprite);
		window.display();
	}
	return 0;
}

int main2()
{
	//fluid parameters
	int xSize = 50;
	int ySize = 50;
	int elementSize = 20;
	float diff = 0.0907;
	float visc = 0.5975;

	//simulation parameters
	int framerate = 60;
	float dt = 1.0 / framerate;
	int numIterations = 5;

	//mouse parameters
	int mouseX = 0;
	int mouseY = 0;
	int lastMouseX = 0;
	int lastMouseY = 0;

	//setting up render window
	sf::RenderWindow window(sf::VideoMode((xSize + 2) * elementSize, (ySize + 2) * elementSize), "Project");
	window.setFramerateLimit(framerate);

	//creating rendering tools
	sf::Image image;
	sf::Texture texture;
	sf::Sprite sprite;
	image.create((xSize + 2) * elementSize, (ySize + 2) * elementSize, sf::Color::Black);

	//now let's create the fluid we want (x, y, diffusion, viscosity, timestep, accuracy/quality)
	auto fluid = (FluidLogic*)(new Fluid(xSize, ySize, diff, visc, dt, numIterations));

	while (window.isOpen())
	{
		//window details
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
				window.close();
		}

		lastMouseX = mouseX;
		lastMouseY = mouseY;
		sf::Vector2i localPosition = sf::Mouse::getPosition(window);
		mouseX = localPosition.x;
		mouseY = localPosition.y;

		if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
			if ((mouseX / elementSize > 0 && mouseX / elementSize < xSize - 1 && mouseY / elementSize  > 0 && mouseY / elementSize < ySize - 1) && (lastMouseX / elementSize > 0 && lastMouseX / elementSize < xSize - 1 && lastMouseY / elementSize  > 0 && lastMouseY / elementSize < ySize - 1)) {
				fluid->AddSourceOfDensity(lastMouseX / elementSize, lastMouseY / elementSize, mouseX / elementSize, mouseY / elementSize);
			}
		}

		if (sf::Mouse::isButtonPressed(sf::Mouse::Right)) {
			fluid->reset();
		}

		fluid->update();
		fluid->render(elementSize, image);
		texture.loadFromImage(image);
		sprite.setTexture(texture, true);
		window.draw(sprite);
		window.display();
	}
	return 0;
}