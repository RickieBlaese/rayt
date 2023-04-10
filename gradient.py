# copied and edited from https://stackoverflow.com/questions/17142331/convert-truetype-glyphs-to-png-image

from PIL import Image, ImageFont, ImageDraw
import sys, math

font_size = 16
font_color = "#FFFFFF"

font = ImageFont.truetype(sys.argv[1], font_size)

desired_characters = open("gradient_source_chars.txt").read()

brightness: dict[str, float] = {}

for character in desired_characters:
    left, top, right, bottom = font.getbbox(character)
    img = Image.new("RGBA", (right - left, bottom - top))
    draw = ImageDraw.Draw(img)
    draw.text((left, top), str(character), font=font, fill=font_color)

    count = 0

    # potential method for comparing to gaussian distribution
    for i in range(img.width):
        for j in range(img.height):
            r, g, b, a = img.getpixel((i, j));
            value = (r + b + g) / (3.0 * 255.0)
            # mean squared error of the color, but just summed up
            count += (math.e ** (- (i / (img.width / 2.0) - 1)**2 - (j / (img.height / 2.0) - 1)**2))**2 - value**2

    # for _, _, _, a in img.getdata():
    #     if a > 0:
    #         count += 1
    brightness[character] = count


brightness = dict(sorted(brightness.items(), key=lambda item: item[1], reverse=True))
file = open("gradient.txt", "w")
for ch in brightness:
    if ch != '\n':
        file.write(ch)

file.close()
