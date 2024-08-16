# Bforartists

Note that the nightlies produced by Github Actions are not meant for productive work. They miss Cuda and Optix. And will most likely not work. The purpose of the nightlies is to  show that the code compiles. Please have a look at the download section at our page. The devbuilds there are working versions.
Goo Engine
==========

Goo Engine is a version of Blender by [DillonGoo Studios](https://www.youtube.com/dillongoo) with some key modifications
to the source code that allow us to focus on NPR and Anime style rendering. Key features include the four custom Shader 
nodes we added to Eevee and Light Groups. You can learn about them on the 
[Professor Goo YouTube channel](https://www.youtube.com/@professorGoo>).

Some have asked why we haven't pushed these changes to Blender's master branch: 
this is because many of these features would not be accepted as they are and may destablize other features like Cycles.
We only use Eevee for NPR so we prioritize specific features which the Blender Institute understandably can't afford to 
make compromises for. This allows us to develop and test quickly. However, we do actively update Goo Engine to keep up 
with the native Blender branch so they can evolve side by side. We also regularly add and improve custom features as we 
use Goo Engine on projects at the studio. [Stay tuned](https://twitter.com/dillongoostudio>) for updates.

Get the pre-built download on the [Patreon page](https://www.patreon.com/dillongoo) to support continued development of 
Goo Engine.


Blender
=======

[![Ubuntu](https://github.com/Bforartists/Bforartists/actions/workflows/linux.yml/badge.svg)](https://github.com/Bforartists/Bforartists/actions/workflows/linux.yml)
[![Mac](https://github.com/Bforartists/Bforartists/actions/workflows/mac.yml/badge.svg)](https://github.com/Bforartists/Bforartists/actions/workflows/mac.yml)
[![Windows](https://github.com/Bforartists/Bforartists/actions/workflows/windows.yml/badge.svg)](https://github.com/Bforartists/Bforartists/actions/workflows/windows.yml)

Bforartists is a fork of the popular 3D software Blender, with the goal to improve the graphical UI, and to increase the usability.

The project page, the binary downloads and further information can be found at the project page: https://www.bforartists.de/

Please use the tracker to add issues, bug reports and requests. The tracker can be found here: https://github.com/Bforartists/Bforartists/issues

![alt text](https://www.bforartists.de/wp-content/uploads/2020/04/modelingbfa2.jpg)

Why choose Bforartists, and not Blender? Because our UI  is streamlined, cleaner, better organized, has colored icons, left aligned text, and much more of all those things that makes a good UI and the life easier. Bforartists is fully compatible with Blender. You can use both asides to make yourself a picture. The files are transferable. So try it, you loose nothing.

The main Differences between Bforartists and Blender are:

- An own keymap, which is reduced to just the necessary hotkeys and a navigation that can be purely done by mouse.
- Cleaned up User Interface. Lots of not necessary double, triple or even more equal menu entries removed.
- Extended User Interface. Lots of formerly hotkey only tools have a menu entry now.
- Rearranged User Interface. Some things are better accessible now, some are not so much in the way anymore.
- Improved defaults.
- Colored and as double as much icons than Blender.
- Configurable UI Elements like the Toolbar, with icon buttons.
- Tabs in the toolshelf.
- Improved Layouts.
- Left aligned checkboxes and text where possible.
- Better Tooltips.
- Better readable standard theme.
- Some neat addons to improve usability, like the reset 3D View addon or the Set Dimensions addon with which you can scale in world coordinates in edit mode.
- And lots more small details like not so much confirm dialogs. Or that we have resurrected the tabs in the tool shelf to allow quicker access to the tools.

A detailed list of the changes can be found in the release notes: www.bforartists.de/wiki/release-notes

But the code is just half of the show. Another important bit are the non code things.

- The target audience for Bforartists are hobbyists and indie developers. Blender tries to target professionals. That’s a completely different audience and development target.
- We have a better manual. The Blender manual is unsearchable in big parts, has an odd structure, and still relies heavily on the Blender keymap to name just a few flaws. This makes it nearly unusable for users. And that’s why we have rewritten it. With the tools in mind, not the hotkeys. Without odd opinions like Angle based is better than LSCM to unwrap. Without general CG Tutorials for a tool, and forgetting to describe how the tool really works and what it does. And a structure that follows the editors and the menus.

Here you can find some examples of the vital differences between Blender and Bforartists: https://www.bforartists.de/the-differences-to-blender/

This youtube video describes the differences by showing them.

[![Youtube Introduction](https://img.youtube.com/vi/xAJQsKRi3sY/0.jpg)](https://www.youtube.com/watch?v=xAJQsKRi3sY)

Here's a Youtube playlist with some quickstart videos:

https://www.youtube.com/watch?v=sZlqqMAGgMs&list=PLB0iqEbIPQTZArhZspyYSJOS_00jURpUB

For those who are willing and interested to help with the manual, we have a Github repository for it too. You can find it here: 

https://github.com/Bforartists/Manual

In case that you want to compile Bforartists by yourself, the advice is the same as for Blender. Here’s a short tutorial that shows the whole process under Linux Ubuntu 17 https://www.bforartists.de/data/tuts/Building_with_Ubuntu17.pdf
And here a youtube video that shows how to compile Bforartists 2 on Ubuntu 20: https://www.youtube.com/watch?v=sUhM2YSw1gQ

Here’s the Blender advice for Linux, which goes a bit more into detail: https://wiki.blender.org/wiki/Building_Blender/Linux
Here’s the advice for Mac: https://wiki.blender.org/index.php/Dev:Doc/Building_Blender/Mac
For Windows have a look here: https://www.bforartists.de/data/tuts/Building_Bforartists_with_Cmake_and_VS_2017_on_Windows.pdf

You might want to compile an older Bforartists version where the newest precompiled libs do not work. Here you can find the precompiled libs for all previous Blender versions:

https://svn.blender.org/svnroot/bf-blender/tags/

Have fun :)

The Bforartists team

What our users says:

”Just want to share some feedback on this version of Blender. So awesome. I finally feel ungimped. Whoever people who made this, you are truly magnificent people. I finally started using it today, and felt let I was using stuff I already knew from Unity and Unreal. It’s a massive relief on the brain not to have to keep translating from one system to the next. I know these software seem like a dime a dozen, but this one serves a very great purpose for people who have accustomed to certain patterns and then having extreme difficulty transferring momentarily to another. You are the greatest.” – Ronin

”superior version of blender” – Keiko Furukawa

”Great work with much enthusiasm!” – Tihomir Dovramadjiev Phd

”la interface es mucho mejor para animar gracias por eso” – Hardy Iglesias Garcia (the interface is much better to animate in, thanks for that)

”If Blender setups give you nightmares, try this.” – George Perkins

”Better than blender. Easy to learn and user friendly” – Emīls Geršinskis-Ješinskis

”5 stars” – Josh Morris

”Great! Keep up the good work!” – Runar Nyegaarden

THANK YOU! I had no idea this existed till yesterday. It’s EXACTLY what I wished for...Blender functionality without the frustratingly convoluted interface. I have stuff to do! I have no problem understanding how the tools work, but Blender gets in my way, and it pisses me off! Geeks step aside...it’s the artists turn to take the controls. Bforartists. Thanks 🙏🍷🎩🎩🎩🎩🎩 - Youtube comment

I’ve been using BforArtists a few months now and absolutely love it. Excellent, excellent job guys! - Youtube comment

This deserves far more attention than it is getting - Youtube comment






