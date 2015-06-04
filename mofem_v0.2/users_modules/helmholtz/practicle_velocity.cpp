


PetscBool practicle_velocity = PETSC_FALSE;
ierr = PetscOptionsGetBool(PETSC_NULL,"-practicle_velocity",&practicle_velocity,NULL); CHKERRQ(ierr);


if(practicle_velocity) {
//Particle Velocity field, do not confused with speed of sound.
ierr = m_field.add_field("reVEL",H1,3); CHKERRQ(ierr);
ierr = m_field.add_field("imVEL",H1,3); CHKERRQ(ierr);
}


if(practicle_velocity) {
//Particle Velocity field, do not confused with speed of sound.
ierr = m_field.add_ents_to_field_by_TETs(root_set,"reVEL"); CHKERRQ(ierr);
ierr = m_field.add_ents_to_field_by_TETs(root_set,"imVEL"); CHKERRQ(ierr);
}



if(practicle_velocity) {
//Particle Velocity field, do not confused with speed of sound.
ierr = m_field.set_field_order(root_set,MBTET,"reVEL",order); CHKERRQ(ierr);
ierr = m_field.set_field_order(root_set,MBTRI,"reVEL",order); CHKERRQ(ierr);
ierr = m_field.set_field_order(root_set,MBEDGE,"reVEL",order); CHKERRQ(ierr);
ierr = m_field.set_field_order(root_set,MBVERTEX,"reVEL",1); CHKERRQ(ierr);
ierr = m_field.set_field_order(root_set,MBTET,"imVEL",order); CHKERRQ(ierr);
ierr = m_field.set_field_order(root_set,MBTRI,"imVEL",order); CHKERRQ(ierr);
ierr = m_field.set_field_order(root_set,MBEDGE,"imVEL",order); CHKERRQ(ierr);
ierr = m_field.set_field_order(root_set,MBVERTEX,"imVEL",1); CHKERRQ(ierr);
}
