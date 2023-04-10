# copied and edited from https://stackoverflow.com/questions/17142331/convert-truetype-glyphs-to-png-image

from PIL import Image, ImageFont, ImageDraw
import sys

font_size = 16
font_color = "#000000"

font = ImageFont.truetype(sys.argv[1], font_size)

desired_characters = open("gradient_source_chars.txt").read()

brightness: dict[str, float] = {}

for character in desired_characters:
    left, top, right, bottom = font.getbbox(character)
    img = Image.new("RGBA", (right - left, bottom - top))
    draw = ImageDraw.Draw(img)
    draw.text((left, top), str(character), font=font, fill=font_color)

    count = 0
    for p in img.getdata():
        if p[3] > 0: # not transparent
            count += 1
    brightness[character] = count


brightness = dict(sorted(brightness.items(), key=lambda item: item[1]))
file = open("gradient.txt", "w")
for ch in brightness:
    if ch != '\n':
        file.write(ch)

file.close()
