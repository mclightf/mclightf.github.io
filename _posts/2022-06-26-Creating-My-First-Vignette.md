
I made my first vignette recently, where I explored this [Pokémon API](https://pokeapi.co/docs/v2). I wanted to chat a little bit about my experience.
For context, feel free to give my vignette a look [here](https://mclightf.github.io/ST558-Project-1/), and the associated github repo [here](https://github.com/mclightf/ST558-Project-1).

In this vignette, I focused on exploring some of the statistics of the moves that Pokémon can learn. I designed a function to query the names, types, powers, accuracies, power points, and damage classes of the moves for any individual Pokémon. In my exploratory data analysis, I focused on the Pokémon Squirtle and its moves. In short, I looked at the distributions of the Power of the moves within the levels of multiple categorical variables and defined some new statistics to help use determine Squirtle's *best* move, which happened to be Water Spout.

**What was the most difficult part of the logic and programming for me?**

The most difficult part of making this vignette for me is definitely figuring out how to get the information that I actually wanted from the Pokémon API. There was no real way around this besides trial and error of reading about the different endpoints, trying different queries, and checking attributes and structures until you find what you need. I also found combining all of this data into a sensible output a bit difficult, since each move for the Pokémon required a separate individual query in order to find the desired statistics, resulting in a pretty substantial loop. There were also some special cases of the Power and Accuracy variables that I had to account for. 

**What would I do differently in approaching a similar project in the future?**

The more work I put into the vignette, the more I realized how much further I kept wanting to take the analysis. Ultimately, I was not able to explore everything that came to mind. I would allow myself more time in the future to truly be able to explore where my analysis may take me. Also, my vignette was mainly a thorough analysis of a specific subset of the Pokémon API, so there were many other endpoints left untouched. I might go a different route next time and explore more endpoints while keeping my analysis more surface-level.

Overall, I had a good time! Expect to see more things like this from me in the future. :)
