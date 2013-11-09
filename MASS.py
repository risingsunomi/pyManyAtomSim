#!/usr/bin/env python

# This a modification of a python program from the tutorial "pygame physics simulation"
# created by Peter Collingridge (@Aleatoric/http://www.petercollingridge.co.uk/) 
# and a many partical simulation Pascal program found in the book 
# "An Introduction to Computer Simulation Methods - Applications to Physical Systems Part 1" 
# by Harvey Gould and Jan Tobochnik

# This simulates interactions of idle gas atoms in a closed vaccume. Temperature,
# pressure and molar concentraction are varibles that can be changed. Area will be
# used instead of volume due to simulation being in 2D 
# Friction, drag, elasticity and gravity is taken out

# risingsunomi

# Friction, drag and elasticity is taken out - dealing on the 10x magnification scale
# see http://chemed.chem.wisc.edu/chempaths/GenChem-Textbook/Kinetic-Theory-of-Gases-Molecular-Speeds-940.html
# for molecular speeds information.

#http://ptp.oxfordjournals.org/content/48/6/2132.full.pdf
#http://gozips.uakron.edu/~mattice/ps674/lj.html

from pygame.locals import *
import pygame
import random
import math

background_colour = (0,0,0)

def addVectors((angle1, length1), (angle2, length2)):
    x  = math.sin(angle1) * length1 + math.sin(angle2) * length2
    y  = math.cos(angle1) * length1 + math.cos(angle2) * length2
    
    angle = 0.5 * math.pi - math.atan2(y, x)
    length  = math.hypot(x, y)

    return (angle, length)

class Environment():
    def __init__(self, W, H, T, P, screen, speed_factor=0.0):
        self.temp = T
        self.pressure = P
        self.width = W
        self.height = H
        self.area = 0.01 * float(self.width * self.height) # for R constant calculations
        self.screen = screen
        self.atoms = []
        self.total_atoms = 0
        self.speed_factor = speed_factor

    def update(self):
        for i, atom in enumerate(self.atoms):
            atom.move()
            atom.bounce()
            for atom2 in self.atoms[i+1:]:
                self.lfeffects(atom, atom2)
                self.collide(atom, atom2)
            atom.display()

    def lfeffects(self, p1, p2):
        # Lennard-Jones potential effects
        dx = p1.x - p2.x
        dy = p1.y - p2.y
        
        dist = math.hypot(dx, dy)
        total_mass = p1.mass + p2.mass

        ljnumber = 4*p1.ljpotential*(((p1.alpha/dist)**12) - ((p1.alpha/dist)**6))
        if ljnumber > 0:
            if p1.charge > p2.charge and p1.mass >= p2.mass:
                p2.angle = p1.angle

            elif p1.charge > p2.charge and p1.mass < p2.mass:
                p1.angle = p2.angle

            elif p1.charge < p2.charge and p1.mass <= p2.mass:
                p1.angle = p2.angle

            elif p1.charge < p2.charge and p1.mass > p2.mass:
                p2.angle = p1.angle

    def collide(self, p1, p2):
        dx = p1.x - p2.x
        dy = p1.y - p2.y
        
        dist = math.hypot(dx, dy)
        if dist < p1.size + p2.size:
            angle = math.atan2(dy, dx) + 0.5 * math.pi
            total_mass = p1.mass + p2.mass

            # inelastic collision and eleastic collision
            # with van der Waals like interactions

            # inelastic collisions due to charge
            if p1.charge != 0 and p2.charge != 0:
                
                if p1.charge > p2.charge and p1.mass >= p2.mass:
                    (p2.angle, p2.speed) = addVectors((p1.angle, (p1.speed*p1.mass)/total_mass), (p1.angle, (p2.speed*p2.mass)/total_mass))

                elif p1.charge > p2.charge and p1.mass < p2.mass:
                    (p1.angle, p1.speed) = addVectors((p2.angle, (p1.speed*p1.mass)/total_mass), (p2.angle, (p2.speed*p2.mass)/total_mass))

                elif p1.charge < p2.charge and p1.mass <= p2.mass:
                    (p2.angle, p2.speed) = addVectors((p2.angle, (p1.speed*p1.mass)/total_mass), (p2.angle, (p2.speed*p2.mass)/total_mass))

                elif p1.charge < p2.charge and p1.mass > p2.mass:
                    (p2.angle, p2.speed) = addVectors((p1.angle, (p1.speed*p1.mass)/total_mass), (p1.angle, (p2.speed*p2.mass)/total_mass))

            # elastic collisions due to no charge
            else:
                (p1.angle, p1.speed) = addVectors((p1.angle, p1.speed*(p1.mass-p2.mass)/total_mass), (angle, 2*p2.speed*p2.mass/total_mass))
                (p2.angle, p2.speed) = addVectors((p2.angle, p2.speed*(p2.mass-p1.mass)/total_mass), (angle+math.pi, 2*p1.speed*p1.mass/total_mass))

            overlap = 0.5*(p1.size + p2.size - dist+1)
            p1.x += math.sin(angle)*overlap
            p1.y -= math.cos(angle)*overlap
            p2.x -= math.sin(angle)*overlap
            p2.y += math.cos(angle)*overlap

    def findatom(self, x, y):
        for p in self.atoms:
            if math.hypot(p.x-x, p.y-y) <= p.size:
                return p
        return None

    #Molecular Speed - speed = sqrt(3RTP/Nm/M)
    #R = 8.314 J/mol K = 8.314 kg m2/s2 mol K
    #R = 0.08205746LatmK-1mol-1
    def mspeed(self, mass):
        # v root-mean square
        if self.temp <= 14.05:
            total_speed = 0
        else:
            total_speed = math.sqrt((3*8.314*(self.temp-14.05))/(float(mass)/1000.00))/self.speed_factor

        return total_speed

class Atom():
    def __init__(self, (x, y), size):
        self.x = x
        self.y = y
        self.size = size
        self.colour = (0, 0, 255)
        self.thickness = 0
        self.speed = 1
        self.angle = math.pi
        self.mass = 1.00794
        self.fpoint = 14.05 # K
        self.charge = random.randint(-8,8)
        self.environment = None
        self.ljpotential = 8.6 #K
        self.alpha = 1 # Angstroms - the distance at which the intermolecular potential between the two particles is zero

    def display(self):
        pygame.draw.circle(self.environment.screen, self.colour, (int(self.x), int(self.y)), self.size, self.thickness)

    def move(self):
        self.x += math.sin(self.angle) * self.environment.mspeed(self.mass)
        self.y -= math.cos(self.angle) * self.environment.mspeed(self.mass)

    def bounce(self):
        if  self.x > self.environment.width - self.size:
            self.x = 2*(self.environment.width - self.size) - self.x
            self.angle = random.uniform(0, math.pi*2) - self.angle

        elif self.x < self.size:
            self.x = 2*self.size - self.x
            self.angle = random.uniform(0, math.pi*2) - self.angle

        if self.y > self.environment.height - self.size:
            self.y = 2*(self.environment.height - self.size) - self.y
            self.angle = math.pi - self.angle

        elif self.y < self.size:
            self.y = 2*self.size - self.y
            self.angle = math.pi - self.angle

def main():
    #flags = FULLSCREEN | DOUBLEBUF
    flags = RESIZABLE | DOUBLEBUF
    screen = pygame.display.set_mode((1120,900), flags)
    #screen = pygame.display.set_mode((0,0), flags)
    screen.set_alpha(None)

    # Set up environment
    cvidinfo = pygame.display.Info()
    cenv = Environment(cvidinfo.current_w, cvidinfo.current_h, 33.15, 0.987, screen, 1000.00) # Pluto 33.15K
    cenv.total_atoms = 250 # 250 and up getting lag - try on better graphics card
    my_atoms = []

    for n in range(cenv.total_atoms):

        # Atom initial configuration
        size = 12 # radius
        x = random.randint(size, cenv.width-size)
        y = random.randint(size, cenv.height-size)

        # Initialize Atom class
        atom = Atom((x, y), size)
        atom.environment = cenv
        cenv.atoms.append(atom)

        atom.colour = (random.uniform(0,184), random.uniform(0,255), random.uniform(0,255))
        atom.speed = cenv.mspeed(atom.mass)
        atom.angle = random.uniform(0, math.pi*2)

        my_atoms.append(atom)

    cenv.atoms = my_atoms
    selected_atom = None
    running = True
    fps_clock = pygame.time.Clock()
    while running:
        fps_tick = fps_clock.tick(40)
        pygame.display.set_caption('Multiple Atoms System Simulator (MASS) v0.1  FPS: %.2f' % fps_clock.get_fps())
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            elif event.type == pygame.MOUSEBUTTONDOWN:
                (mouseX, mouseY) = pygame.mouse.get_pos()
                selected_atom = cenv.findatom(mouseX, mouseY)
                print mouseX, mouseY, selected_atom
            elif event.type == pygame.MOUSEBUTTONUP:
                selected_atom = None
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_ESCAPE:
                    running = False
            elif event.type == VIDEORESIZE: # simulate pressure change by changing volume
                screen = pygame.display.set_mode(event.dict['size'], RESIZABLE | DOUBLEBUF)
                cenv.screen = screen

                cvidinfo = pygame.display.Info()
                cenv.width = cvidinfo.current_w
                cenv.height = cvidinfo.current_h 
                
                for atom in cenv.atoms:
                    atom.environment = cenv
                cenv.update()

        if selected_atom:
            (mouseX, mouseY) = pygame.mouse.get_pos()
            dx = mouseX - selected_atom.x
            dy = mouseY - selected_atom.y
            selected_atom.angle = 0.5*math.pi + math.atan2(dy, dx)
            selected_atom.speed = (math.hypot(dx, dy) * 0.1)*cenv.speed_factor

        screen.fill(background_colour)
        cenv.update()
        
        pygame.display.flip()

if __name__ == "__main__":
    main()