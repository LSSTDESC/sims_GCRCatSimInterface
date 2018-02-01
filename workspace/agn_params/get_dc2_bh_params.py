import GCRCatalogs

if __name__ == "__main__":

    cat = GCRCatalogs.load_catalog('proto-dc2_v2.1.2')

    qties = cat.get_quantities(['redshift_true', 'blackHoleMass',
                                'blackHoleAccretionRate'])

    with open('data/proto_dc2_bh_params.txt', 'w') as out_file:
        for i_obj in range(len(qties['redshift_true'])):
            if qties['blackHoleMass'][i_obj]>0.0:
                if qties['blackHoleAccretionRate'][i_obj]>0.0:
                    out_file.write('%e %e %e\n' %
                                   (qties['blackHoleMass'][i_obj],
                                    qties['blackHoleAccretionRate'][i_obj],
                                    qties['redshift_true'][i_obj]))
