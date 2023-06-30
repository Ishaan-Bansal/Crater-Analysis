import csv
import os

def convert_to_csv(folder_path, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Get a list of all files in the folder
        files = os.listdir(folder_path)
        
        # Loop through each file
        for file_name in files:
            if file_name.endswith('.txt'):
                file_path = os.path.join(folder_path, file_name)
                
                # Read the file contents
                with open(file_path, 'r') as file:
                    content = file.read()
                    
                    # Split the content by commas
                    parts = content.split(',')
                    
                    # Remove the first part (before the first colon)
                    cleaned_parts = [file_name]
                    for part in parts:
                        split_part = part.split(':', 1)
                        if len(split_part) > 1:
                            cleaned_parts.append(split_part[1].strip())
                    
                    # Write the cleaned parts as a row in the CSV file
                    writer.writerow(cleaned_parts)
    
    print('Conversion completed successfully.')

# Provide the folder path containing the text files
folder_path = 'Lab Craters\Combined Results\TXTs'

# Provide the output CSV file path
output_file = 'Lab Craters\Combined Results/analysis.csv'

# Call the function to convert the text files to CSV
convert_to_csv(folder_path, output_file)
