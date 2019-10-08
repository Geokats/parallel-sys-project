import csv
import sys

if __name__ == '__main__':
    in_file = sys.argv[1]
    spd_file = sys.argv[2]
    eff_file = sys.argv[3]
    print "Reading from file %s and writing to files %s and %s" % (in_file, spd_file, eff_file)

    #Get rows
    with open(in_file, mode='r') as f:
        csv_r = csv.reader(f)
        rows = []
        for row in csv_r:
            rows.append(row)

    #Calculate speed-up
    with open(spd_file, mode='w') as spd_f, open(eff_file, mode='w') as eff_f:
        csv_spd = csv.writer(spd_f, delimiter=',')
        csv_spd.writerow(rows[0])

        csv_eff = csv.writer(eff_f, delimiter=',')
        csv_eff.writerow(rows[0])

        for i in range(1, len(rows)):
            spd_row = [rows[i][0], rows[i][1]]

            eff_row = [rows[i][0], rows[i][1]]
            for j in range(2, len(rows[i])):
                try:
                    spd_val = round(float(rows[1][j])/float(rows[i][j]), 6)

                    eff_val = round(float(rows[1][j])/(float(rows[i][j]) * float(rows[i][0])), 6)
                except Exception as e:
                    spd_val = 0

                    eff_val = 0
                spd_val = "{0:0.2f}".format(spd_val)
                spd_row.append(str(spd_val))

                eff_val = "{0:0.2f}".format(eff_val)
                eff_row.append(str(eff_val))

            csv_spd.writerow(spd_row)

            csv_eff.writerow(eff_row)
