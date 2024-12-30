import requests
from bs4 import BeautifulSoup
from collections import defaultdict
import pandas as pd
import re
import time
import os

#%% Download Asteroid Model from 3D Asteroids Website

def download_obj_file(asteroid_id, Asteroid_Name, log_callback):
    first_Letter = Asteroid_Name[0].lower()
    base_url = f"https://3d-asteroids.space/data/asteroids/models/{first_Letter}/{asteroid_id}_{Asteroid_Name}"
    base_dir = 'Asteroids/'
    
    # First, try the untagged file
    url = f"{base_url}.obj"
    # log_callback(f"Trying URL: {url}")
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            directory = os.path.join(base_dir, Asteroid_Name)
            if not os.path.exists(directory):
                os.makedirs(directory)
            file_path = os.path.join(directory, f"{Asteroid_Name}.obj")
            with open(file_path, 'wb') as file:
                file.write(response.content)
            log_callback(f"Downloaded {asteroid_id} model to {file_path}")
        elif response.status_code != 404:
            log_callback(f"Failed to download {asteroid_id} model. Status code: {response.status_code}")
    except requests.exceptions.RequestException as e:
        log_callback(f"Request failed: {e}")
    
    # Try tagged files starting from 1000
    for i in range(10, 5000):  # Start at 1000 and go up to 4999
        url = f"{base_url}_{i}.obj"
        
        # log_callback(f"Trying URL: {url}")
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                directory = os.path.join(base_dir, f"{Asteroid_Name}_{i}")
                if not os.path.exists(directory):
                    os.makedirs(directory)
                file_path = os.path.join(directory, f"{Asteroid_Name}_{i}.obj")
                with open(file_path, 'wb') as file:
                    file.write(response.content)
                log_callback(f"Downloaded {asteroid_id} model to {file_path}")

        except requests.exceptions.RequestException as e:
            log_callback(f"Request failed: {e}")
        
        time.sleep(0.5)  # Add a small delay between requests

    log_callback(f"Finished checking all tags for {asteroid_id}.")


# # Example usage
# asteroid_id = '33181'  # Replace with the actual asteroid ID
# output_path = 'Aalokpatwa.obj'
# Asteroid = 'Aalokpatwa'

# download_obj_file(asteroid_id, output_path, Asteroid)

# Example usage
# asteroid_id = '99942'  # Replace with the actual asteroid ID
# output_path = 'Apophis.obj'
# Asteroid = 'Apophis'

# download_obj_file(asteroid_id, output_path, Asteroid)

#%% List ASteroids From 3D Asteroids Website
def list_asteroids():
    url = 'https://3d-asteroids.space/asteroids/'
    response = requests.get(url)
    
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        # Find all divs with class 'blk'
        asteroid_divs = soup.find_all('div', class_='blk')
        
    asteroids_dict = defaultdict(list)
    named_started = False

    for div in asteroid_divs:
        asteroid_links = div.find_all('a', href=True)
        asteroids = [link.text.strip() for link in asteroid_links if link.text.strip() and '▲Top' not in link.text]
        for asteroid in asteroids:
            # Extract number in parentheses
            match = re.match(r'^\((\d+)\)\s*(.*)', asteroid)
            if match:
                number = match.group(1)
                stripped_name = match.group(2)
            else:
                number = None
                stripped_name = asteroid

            first_letter = stripped_name[0].upper()
            if first_letter.isalpha():
                named_started = True
                asteroids_dict[first_letter].append((stripped_name, number))
            else:
                if not named_started:
                    asteroids_dict['Unnamed'].append((stripped_name, number))
                else:
                    asteroids_dict[first_letter].append((stripped_name, number))

    # Convert the dictionary to a DataFrame
    asteroid_list = []
    for letter, names in asteroids_dict.items():
        for name, number in names:
            asteroid_list.append({'Category': letter, 'Name': name, 'ID': number})
    
    df = pd.DataFrame(asteroid_list)
    return df

# # Example usage

# df_asteroids = list_asteroids(url)

# # Display DataFrame
# df_asteroids.head()


