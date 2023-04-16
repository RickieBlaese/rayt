# copied and edited from https://stackoverflow.com/questions/17142331/convert-truetype-glyphs-to-png-image

from PIL import Image, ImageFont, ImageDraw
import sys, math

font_size = 50
font_color = "#FFFFFF"

font = ImageFont.truetype(sys.argv[1], font_size)

desired_characters = open("gradient_source_chars.txt").read()

brightness: dict[str, float] = {}

def stddev(l: list[float]) -> float:
    mean = 0
    for i in l: mean += i
    mean /= len(l)
    var = 0
    for i in l: var += (i - mean) ** 2
    var /= len(l)
    return math.sqrt(var)

def partition(a: int, b: int, count: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    if count <= 1:
        out.append((a, b))
        return out
    per: float = (b - a) / count
    i: int = 0
    while i < count - 1:
        out.append((int(i * per + a), int((i + 1) * per + a)))
        i += 1
    out.append((int(i * per + a), b))
    return out

for character in desired_characters:
    left, top, right, bottom = font.getbbox(character)
    img = Image.new("RGBA", (right, bottom))
    draw = ImageDraw.Draw(img)
    draw.text((0, 0), str(character), font=font, fill=font_color)

    count = 0

    for r, g, b, a in img.getdata():
        count += a

    # for i in range(img.width):
    #     for j in range(img.height):
    #         r, g, b, a = img.getpixel((i, j));
    #         value = (r + b + g) / (3.0 * 255.0)
    #         # mean squared error of the color, but just summed up
    #         count += (math.e ** (- (i / (img.width / 2.0) - 1)**2 - (j / (img.height / 2.0) - 1)**2))**2 - value**2


    # weird smoothing method
    # l: list[float] = []
    # xcoords = partition(0, img.width, width_block_count)
    # ycoords = partition(0, img.height, height_block_count)
    # for ax, bx in xcoords:
    #     for ay, by in ycoords:
    #         mean = 0
    #         for i in range(ax, bx):
    #             for j in range(ay, by):
    #                 r, b, g, _ = img.getpixel((i, j))
    #                 value = (r + b + g) / 3.0
    #                 mean += value
    #                 count += value
    #         tw = bx - ax
    #         th = by - ay
    #         if abs(tw) < 0.01 or abs(th) < 0.01:
    #             print(xcoords, ycoords, ax, bx, ay, by, mean)
    #         l.append(mean / (tw * th))
    # sdev = stddev(l)

    brightness[character] = count


brightness = dict(sorted(brightness.items(), key=lambda item: item[1]))
file = open("gradient.txt", "w")
for ch in brightness:
    if ch != '\n':
        file.write(ch)

file.close()
